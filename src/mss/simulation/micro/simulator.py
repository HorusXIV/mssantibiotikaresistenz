"""
MicroSimulator - Interface for batch processing of patient episodes.

Designed for parallel execution across thousands of patients.
"""

from __future__ import annotations

import multiprocessing as mp
from concurrent.futures import Executor, ProcessPoolExecutor, ThreadPoolExecutor
from math import ceil
from typing import Any, Dict, Iterable, List, Optional

import numpy as np

from .config import SimulationConfig
from .engine import (
    StrainPopulation,
    build_founder_pool,
    population_to_response,
    simulate_day,
)
from .state import EpisodeState


def _create_initial_population(
    *,
    episode_id: str,
    dominant_genotype: str,
    dominant_strain_name: str | None,
    resistant_fraction: float,
    config: SimulationConfig,
    seed: int | None,
    seed_genome: np.ndarray | None = None,
) -> StrainPopulation:
    """Create the starting strain population for a micro episode.

    Args:
        episode_id: Episode identifier used as the strain namespace.
        dominant_genotype: Target dominant genotype for initialization.
        dominant_strain_name: Optional preferred name for the dominant strain.
        resistant_fraction: Initial resistant population fraction.
        config: Complete micro simulation configuration.
        seed: Optional random seed for initialization.
        seed_genome: Optional transmitted genome to seed the population.

    Returns:
        Initialized strain population.
    """
    rng = np.random.default_rng(seed)
    founder_pool = build_founder_pool(config) if seed_genome is None else []
    return StrainPopulation.create_initial(
        resistant_fraction=resistant_fraction,
        dominant_genotype=dominant_genotype,
        rng=rng,
        dominant_strain_name=dominant_strain_name,
        founder_pool=founder_pool,
        strain_namespace=episode_id or "episode",
        seed_genome=seed_genome,
    )


def _process_single_request(args: tuple) -> Dict[str, Any]:
    """Process one micro request in a worker-compatible form.

    Args:
        args: Tuple containing the request dictionary, optional episode-state payload,
            and serialized micro configuration.

    Returns:
        Response dictionary compatible with ``Patient.apply_micro_response``. The private
        ``_state`` key contains the updated episode payload for the owning simulator.

    Raises:
        ValueError: If the worker configuration payload is missing.
    """
    request, state_dict, config_dict = args

    if not config_dict:
        raise ValueError("Worker micro config payload is required.")
    config = SimulationConfig(**config_dict)

    # Reconstruct or create population
    if state_dict is not None:
        state = EpisodeState.from_payload(state_dict)
        population = state.population
    else:
        # New episode - create initial population
        initial = request.get("initial_state", {})
        resistant_fraction = initial.get("resistant_fraction", 0.0)
        dominant_genotype = initial.get("dominant_genotype", "S")
        raw_seed_genome = initial.get("seed_genome")
        seed_genome = (
            np.array(raw_seed_genome, dtype=np.float32) if raw_seed_genome is not None else None
        )
        population = _create_initial_population(
            episode_id=str(request.get("episode_id", "")),
            dominant_genotype=dominant_genotype,
            dominant_strain_name=initial.get("dominant_strain_name"),
            resistant_fraction=resistant_fraction,
            config=config,
            seed=request.get("seed", None),
            seed_genome=seed_genome,
        )

    # Extract parameters from request
    abx = request.get("abx", {})
    abx_class = abx.get("class", "none") if abx.get("on", False) else "none"
    dose_level = abx.get("dose_level", "std")
    adherence = request.get("adherence", 1.0)

    host = request.get("host", {})
    immune_strength = host.get("immune_strength", 1.0)

    # Run simulation
    final_pop, history = simulate_day(
        population=population,
        abx_class=abx_class,
        dose_level=dose_level,
        adherence=adherence,
        immune_strength=immune_strength,
        config=config,
        seed=request.get("seed", None),
    )

    # Build response
    response = population_to_response(final_pop, immune_strength=immune_strength, config=config)

    # Add episode tracking
    response["episode_id"] = request.get("episode_id")
    response["patient_id"] = request.get("patient_id")
    response["t_day"] = request.get("t_day", 0)

    # Include updated state for persistence
    response["_state"] = EpisodeState(
        episode_id=request.get("episode_id", ""),
        patient_id=request.get("patient_id", ""),
        population=final_pop,
        day=request.get("t_day", 0) + 1,
    ).to_payload()

    return response


def _process_request_chunk(chunk: List[tuple]) -> List[Dict[str, Any]]:
    """Process one chunk of worker request tuples.

    Args:
        chunk: Request tuples accepted by ``_process_single_request``.

    Returns:
        Response dictionaries in chunk order.
    """
    return [_process_single_request(args) for args in chunk]


class MicroSimulator:
    """Stateful batch interface for within-host micro simulation.

    The simulator owns active episode states, serializes requests for worker execution,
    and stores updated state payloads after each processed request.
    """

    def __init__(self, config: SimulationConfig, n_workers: int = None):
        """Initialize the micro simulator.

        Args:
            config: Complete micro simulation configuration.
            n_workers: Number of parallel workers. If omitted, CPU count is used.
        """
        self.config = config
        self.n_workers = n_workers or mp.cpu_count()
        self._episode_states: Dict[str, EpisodeState] = {}
        self._executor: Executor | None = None

    def _config_payload(self) -> Dict[str, Any]:
        """Serialize the active simulation config for worker dispatch.

        Returns:
            Dictionary containing all ``SimulationConfig`` fields.
        """
        return {
            "steps_per_day": self.config.steps_per_day,
            "max_strains": self.config.max_strains,
            "carrying_capacity": self.config.carrying_capacity,
            "min_population": self.config.min_population,
            "clearance_threshold": self.config.clearance_threshold,
            "base_mutation_rate": self.config.base_mutation_rate,
            "mutation_std": self.config.mutation_std,
            "stress_mutation_boost": self.config.stress_mutation_boost,
            "mutant_transfer_fraction": self.config.mutant_transfer_fraction,
            "base_hgt_rate": self.config.base_hgt_rate,
            "hgt_gene_transfer_prob": self.config.hgt_gene_transfer_prob,
            "selection_strength": self.config.selection_strength,
            "abx_selection_pressure_multiplier": self.config.abx_selection_pressure_multiplier,
            "antibiotic_kill_scale": self.config.antibiotic_kill_scale,
            "growth_rate_per_step": self.config.growth_rate_per_step,
            "death_rate_per_step": self.config.death_rate_per_step,
            "death_fitness_floor": self.config.death_fitness_floor,
            "strain_prune_threshold": self.config.strain_prune_threshold,
            "base_damage_per_step": self.config.base_damage_per_step,
            "replication_damage_factor": self.config.replication_damage_factor,
            "stress_damage_factor": self.config.stress_damage_factor,
            "repair_rate_per_step": self.config.repair_rate_per_step,
            "age_mortality_scale": self.config.age_mortality_scale,
            "damage_mortality_scale": self.config.damage_mortality_scale,
            "lifecycle_half_life_steps": self.config.lifecycle_half_life_steps,
            "max_damage_load": self.config.max_damage_load,
            "dormancy_growth_penalty": self.config.dormancy_growth_penalty,
            "synergy_repair_dormancy_bonus": self.config.synergy_repair_dormancy_bonus,
            "synergy_stress_tolerance_bonus": self.config.synergy_stress_tolerance_bonus,
            "stochastic_threshold": self.config.stochastic_threshold,
            "stochastic_noise_scale": self.config.stochastic_noise_scale,
            "founder_pool_size": self.config.founder_pool_size,
            "founder_pool_seed": self.config.founder_pool_seed,
            "founder_pool_gene_noise_std": self.config.founder_pool_gene_noise_std,
            "gene_presence_threshold": self.config.gene_presence_threshold,
        }

    def _get_executor(self) -> Executor:
        """Lazily create the worker executor.

        Returns:
            Process pool when available, otherwise a thread pool fallback.
        """
        if self._executor is None:
            try:
                self._executor = ProcessPoolExecutor(max_workers=self.n_workers)
            except (OSError, PermissionError):
                # Some restricted runtimes disallow process synchronization primitives.
                self._executor = ThreadPoolExecutor(max_workers=self.n_workers)
        return self._executor

    def _build_args_list(self, requests: List[Dict[str, Any]]) -> List[tuple]:
        """Attach current episode state and config to each request.

        Args:
            requests: Micro request dictionaries from macro patients.

        Returns:
            Worker argument tuples in request order.

        Raises:
            ValueError: If the same episode ID appears more than once in one batch.
        """
        config_payload = self._config_payload()
        args_list = []

        seen_episode_ids = set()
        for req in requests:
            episode_id = req.get("episode_id")
            if episode_id:
                if episode_id in seen_episode_ids:
                    raise ValueError(
                        f"Duplicate episode_id '{episode_id}' in one batch is not supported."
                    )
                seen_episode_ids.add(episode_id)

            state_payload = None
            if episode_id and episode_id in self._episode_states:
                state_payload = self._episode_states[episode_id].to_payload()

            args_list.append((req, state_payload, config_payload))

        return args_list

    def _chunk_args(self, args_list: List[tuple]) -> List[List[tuple]]:
        """Group worker arguments into chunks.

        Args:
            args_list: Worker argument tuples.

        Returns:
            List of argument chunks sized for the configured worker count.
        """
        if len(args_list) <= 1:
            return [args_list]

        target_chunks = min(self.n_workers, len(args_list))
        chunk_size = max(1, ceil(len(args_list) / target_chunks))
        return [args_list[i : i + chunk_size] for i in range(0, len(args_list), chunk_size)]

    def _store_responses(self, responses: Iterable[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Update episode state from responses and remove private payloads.

        Args:
            responses: Response dictionaries produced by worker processing.

        Returns:
            Cleaned responses without the private ``_state`` key.
        """
        cleaned = []
        for resp in responses:
            episode_id = resp.get("episode_id")
            if episode_id and "_state" in resp:
                self._episode_states[episode_id] = EpisodeState.from_payload(resp["_state"])
            resp.pop("_state", None)
            cleaned.append(resp)
        return cleaned

    def process_request(self, request: Dict[str, Any]) -> Dict[str, Any]:
        """Process a single micro request.

        Args:
            request: Dictionary produced by ``Patient.make_micro_request``.

        Returns:
            Response dictionary for ``Patient.apply_micro_response``.
        """
        response = _process_single_request(self._build_args_list([request])[0])
        return self._store_responses([response])[0]

    def process_batch(
        self, requests: List[Dict[str, Any]], parallel: bool = True
    ) -> List[Dict[str, Any]]:
        """Process multiple micro requests.

        Args:
            requests: Request dictionaries from macro patients.
            parallel: Whether to use the worker executor when more than one request is present.

        Returns:
            Response dictionaries in the same order as ``requests``.
        """
        if not requests:
            return []

        args_list = self._build_args_list(requests)

        # Process
        if parallel and len(requests) > 1 and self.n_workers > 1:
            chunks = self._chunk_args(args_list)
            chunked_responses = list(self._get_executor().map(_process_request_chunk, chunks))
            responses = [resp for chunk in chunked_responses for resp in chunk]
        else:
            responses = [_process_single_request(args) for args in args_list]

        return self._store_responses(responses)

    def clear_episode(self, episode_id: str) -> None:
        """Remove persisted state for an episode.

        Args:
            episode_id: Episode identifier to remove.
        """
        self._episode_states.pop(episode_id, None)

    def initialize_episode(
        self,
        *,
        episode_id: str,
        patient_id: str,
        resistant_fraction: float = 0.0,
        dominant_genotype: str = "S",
        dominant_strain_name: str | None = None,
        seed_genome: np.ndarray | None = None,
        day: int = 0,
        seed: int | None = None,
    ) -> EpisodeState:
        """Create and persist an initial episode state without simulating a day.

        Args:
            episode_id: Episode identifier to initialize.
            patient_id: Patient identifier associated with the episode.
            resistant_fraction: Initial resistant population fraction.
            dominant_genotype: Target dominant genotype for initialization.
            dominant_strain_name: Optional preferred dominant strain name.
            seed_genome: Optional transmitted genome to seed the episode.
            day: Initial macro day stored in the episode state.
            seed: Optional random seed for initialization.

        Returns:
            Existing or newly created episode state.
        """
        existing = self._episode_states.get(episode_id)
        if existing is not None:
            return existing

        population = _create_initial_population(
            episode_id=episode_id,
            dominant_genotype=dominant_genotype,
            dominant_strain_name=dominant_strain_name,
            resistant_fraction=resistant_fraction,
            config=self.config,
            seed=seed,
            seed_genome=seed_genome,
        )
        state = EpisodeState(
            episode_id=episode_id,
            patient_id=patient_id,
            population=population,
            day=day,
        )
        self._episode_states[episode_id] = state
        return state

    def get_episode_state(self, episode_id: str) -> Optional[EpisodeState]:
        """Return the current state for an episode.

        Args:
            episode_id: Episode identifier to look up.

        Returns:
            Episode state if active, otherwise ``None``.
        """
        return self._episode_states.get(episode_id)

    def get_active_episodes(self) -> List[str]:
        """Return active episode identifiers.

        Returns:
            List of episode IDs currently tracked by this simulator.
        """
        return list(self._episode_states.keys())

    def close(self) -> None:
        """Shut down the background worker executor if it exists.

        Returns:
            None.
        """
        if self._executor is not None:
            self._executor.shutdown(wait=True)
            self._executor = None

    def __del__(self):
        """Best-effort cleanup for worker resources during garbage collection."""
        try:
            self.close()
        except Exception:
            pass


# Convenience function for standalone use
def run_micro_simulation(
    patient_request: Dict[str, Any], config: SimulationConfig
) -> Dict[str, Any]:
    """Run one stateless micro simulation from a patient request.

    Args:
        patient_request: Dictionary produced by ``Patient.make_micro_request``.
        config: Complete micro simulation configuration.

    Returns:
        Response dictionary for ``Patient.apply_micro_response``. No episode state is
        retained between calls.
    """
    # Create initial population
    initial = patient_request.get("initial_state", {})
    resistant_fraction = initial.get("resistant_fraction", 0.0)
    dominant_genotype = initial.get("dominant_genotype", "S")
    raw_seed_genome = initial.get("seed_genome")
    seed_genome = (
        np.array(raw_seed_genome, dtype=np.float32) if raw_seed_genome is not None else None
    )

    population = _create_initial_population(
        episode_id=str(patient_request.get("episode_id", "")),
        dominant_genotype=dominant_genotype,
        dominant_strain_name=initial.get("dominant_strain_name"),
        resistant_fraction=resistant_fraction,
        config=config,
        seed=patient_request.get("seed", None),
        seed_genome=seed_genome,
    )

    # Extract parameters
    abx = patient_request.get("abx", {})
    abx_class = abx.get("class", "none") if abx.get("on", False) else "none"
    dose_level = abx.get("dose_level", "std")
    adherence = patient_request.get("adherence", 1.0)

    host = patient_request.get("host", {})
    immune_strength = host.get("immune_strength", 1.0)

    # Run simulation
    final_pop, _ = simulate_day(
        population=population,
        abx_class=abx_class,
        dose_level=dose_level,
        adherence=adherence,
        immune_strength=immune_strength,
        config=config,
        seed=patient_request.get("seed", None),
    )

    # Build response
    response = population_to_response(final_pop, immune_strength=immune_strength, config=config)

    response["episode_id"] = patient_request.get("episode_id")
    response["patient_id"] = patient_request.get("patient_id")

    return response
