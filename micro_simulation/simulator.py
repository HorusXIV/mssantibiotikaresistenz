"""
MicroSimulator - Interface for batch processing of patient episodes.

Designed for parallel execution across thousands of patients.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Dict, Any, List, Optional
from concurrent.futures import ThreadPoolExecutor
import multiprocessing as mp

from .simulation import (
    SimulationConfig,
    StrainPopulation,
    simulate_day,
    population_to_response,
)
from .genome import NUM_GENES


@dataclass
class EpisodeState:
    """Persistent state for a patient episode (carried between days)."""

    episode_id: str
    patient_id: str
    population: StrainPopulation
    day: int = 0

    def to_dict(self) -> Dict[str, Any]:
        """Serialize for storage."""
        return {
            "episode_id": self.episode_id,
            "patient_id": self.patient_id,
            "day": self.day,
            "genomes": self.population.genomes.tolist(),
            "populations": self.population.populations.tolist(),
            "lineage_ages": self.population.lineage_ages.tolist(),
            "damage_loads": self.population.damage_loads.tolist(),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> EpisodeState:
        """Deserialize from storage."""
        genomes = np.array(data["genomes"], dtype=np.float32)
        if genomes.size == 0:
            genomes = genomes.reshape(0, NUM_GENES)
        populations = np.array(data["populations"], dtype=np.float64)
        lineage_ages = np.array(data.get("lineage_ages", []), dtype=np.float64)
        damage_loads = np.array(data.get("damage_loads", []), dtype=np.float64)
        return cls(
            episode_id=data["episode_id"],
            patient_id=data["patient_id"],
            population=StrainPopulation(
                genomes=genomes,
                populations=populations,
                lineage_ages=lineage_ages if len(lineage_ages) else None,
                damage_loads=damage_loads if len(damage_loads) else None,
            ),
            day=data.get("day", 0),
        )


def _process_single_request(args: tuple) -> Dict[str, Any]:
    """
    Process a single micro request (for parallel execution).

    Args:
        args: Tuple of (request_dict, episode_state_dict, config_dict)

    Returns:
        Response dict compatible with Patient.apply_micro_response()
    """
    request, state_dict, config_dict = args

    config = SimulationConfig(**config_dict) if config_dict else SimulationConfig()

    # Reconstruct or create population
    if state_dict is not None:
        state = EpisodeState.from_dict(state_dict)
        population = state.population
    else:
        # New episode - create initial population
        initial = request.get("initial_state", {})
        resistant_fraction = initial.get("resistant_fraction", 0.0)

        rng = np.random.default_rng(request.get("seed", None))
        population = StrainPopulation.create_initial(resistant_fraction=resistant_fraction, rng=rng)

    # Extract parameters from request
    abx = request.get("abx", {})
    abx_class = abx.get("class", "none") if abx.get("on", False) else "none"
    dose_level = abx.get("dose_level", "std")
    adherence = request.get("adherence", 1.0)

    host = request.get("host", {})
    immune_strength = host.get("immune_strength", 1.0)
    immune_status = host.get("immune_status", "normal")

    # Run simulation
    final_pop, history = simulate_day(
        population=population,
        abx_class=abx_class,
        dose_level=dose_level,
        adherence=adherence,
        immune_strength=immune_strength,
        immune_status=immune_status,
        config=config,
        seed=request.get("seed", None),
    )

    # Build response
    response = population_to_response(
        final_pop, immune_strength=immune_strength, immune_status=immune_status, config=config
    )

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
    ).to_dict()

    return response


class MicroSimulator:
    """
    High-level interface for micro-simulation.

    Manages episode states and supports batch processing.
    """

    def __init__(self, config: SimulationConfig = None, n_workers: int = None):
        """
        Initialize simulator.

        Args:
            config: Simulation configuration
            n_workers: Number of parallel workers (None = CPU count)
        """
        self.config = config or SimulationConfig()
        self.n_workers = n_workers or mp.cpu_count()
        self._episode_states: Dict[str, EpisodeState] = {}

    def process_request(self, request: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process a single micro request.

        Args:
            request: Dict from Patient.make_micro_request()

        Returns:
            Response dict for Patient.apply_micro_response()
        """
        episode_id = request.get("episode_id")

        # Get existing state if available
        state_dict = None
        if episode_id and episode_id in self._episode_states:
            state_dict = self._episode_states[episode_id].to_dict()

        config_dict = {
            "steps_per_day": self.config.steps_per_day,
            "max_strains": self.config.max_strains,
            "carrying_capacity": self.config.carrying_capacity,
            "min_population": self.config.min_population,
            "clearance_threshold": self.config.clearance_threshold,
            "base_mutation_rate": self.config.base_mutation_rate,
            "mutation_std": self.config.mutation_std,
            "stress_mutation_boost": self.config.stress_mutation_boost,
            "base_hgt_rate": self.config.base_hgt_rate,
            "hgt_gene_transfer_prob": self.config.hgt_gene_transfer_prob,
            "selection_strength": self.config.selection_strength,
            "growth_rate_per_step": self.config.growth_rate_per_step,
            "death_rate_per_step": self.config.death_rate_per_step,
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
        }

        response = _process_single_request((request, state_dict, config_dict))

        # Update stored state
        if episode_id and "_state" in response:
            self._episode_states[episode_id] = EpisodeState.from_dict(response["_state"])

        # Remove internal state from response
        response.pop("_state", None)

        return response

    def process_batch(
        self, requests: List[Dict[str, Any]], parallel: bool = True
    ) -> List[Dict[str, Any]]:
        """
        Process multiple requests, optionally in parallel.

        Args:
            requests: List of request dicts
            parallel: Whether to use parallel processing

        Returns:
            List of response dicts (same order as requests)
        """
        if not requests:
            return []

        # Prepare arguments
        config_dict = {
            "steps_per_day": self.config.steps_per_day,
            "max_strains": self.config.max_strains,
            "carrying_capacity": self.config.carrying_capacity,
            "min_population": self.config.min_population,
            "clearance_threshold": self.config.clearance_threshold,
            "base_mutation_rate": self.config.base_mutation_rate,
            "mutation_std": self.config.mutation_std,
            "stress_mutation_boost": self.config.stress_mutation_boost,
            "base_hgt_rate": self.config.base_hgt_rate,
            "hgt_gene_transfer_prob": self.config.hgt_gene_transfer_prob,
            "selection_strength": self.config.selection_strength,
            "growth_rate_per_step": self.config.growth_rate_per_step,
            "death_rate_per_step": self.config.death_rate_per_step,
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
        }

        args_list = []
        for req in requests:
            episode_id = req.get("episode_id")
            state_dict = None
            if episode_id and episode_id in self._episode_states:
                state_dict = self._episode_states[episode_id].to_dict()
            args_list.append((req, state_dict, config_dict))

        # Process
        if parallel and len(requests) > 1:
            # Use ThreadPoolExecutor for GIL-released NumPy operations
            # For true parallelism with heavy computation, use ProcessPoolExecutor
            # but that requires serialization overhead
            with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
                responses = list(executor.map(_process_single_request, args_list))
        else:
            responses = [_process_single_request(args) for args in args_list]

        # Update states
        for resp in responses:
            episode_id = resp.get("episode_id")
            if episode_id and "_state" in resp:
                self._episode_states[episode_id] = EpisodeState.from_dict(resp["_state"])
            resp.pop("_state", None)

        return responses

    def clear_episode(self, episode_id: str) -> None:
        """Remove episode state (e.g., when patient clears infection)."""
        self._episode_states.pop(episode_id, None)

    def get_episode_state(self, episode_id: str) -> Optional[EpisodeState]:
        """Get current state for an episode."""
        return self._episode_states.get(episode_id)

    def get_active_episodes(self) -> List[str]:
        """Get list of active episode IDs."""
        return list(self._episode_states.keys())


# Convenience function for standalone use
def run_micro_simulation(
    patient_request: Dict[str, Any], config: SimulationConfig = None
) -> Dict[str, Any]:
    """
    Run a single micro simulation from a patient request.

    Stateless - does not persist episode state.

    Args:
        patient_request: Dict from Patient.make_micro_request()
        config: Optional simulation config

    Returns:
        Response dict for Patient.apply_micro_response()
    """
    config = config or SimulationConfig()

    # Create initial population
    initial = patient_request.get("initial_state", {})
    resistant_fraction = initial.get("resistant_fraction", 0.0)

    rng = np.random.default_rng(patient_request.get("seed", None))
    population = StrainPopulation.create_initial(resistant_fraction=resistant_fraction, rng=rng)

    # Extract parameters
    abx = patient_request.get("abx", {})
    abx_class = abx.get("class", "none") if abx.get("on", False) else "none"
    dose_level = abx.get("dose_level", "std")
    adherence = patient_request.get("adherence", 1.0)

    host = patient_request.get("host", {})
    immune_strength = host.get("immune_strength", 1.0)
    immune_status = host.get("immune_status", "normal")

    # Run simulation
    final_pop, _ = simulate_day(
        population=population,
        abx_class=abx_class,
        dose_level=dose_level,
        adherence=adherence,
        immune_strength=immune_strength,
        immune_status=immune_status,
        config=config,
        seed=patient_request.get("seed", None),
    )

    # Build response
    response = population_to_response(
        final_pop, immune_strength=immune_strength, immune_status=immune_status, config=config
    )

    response["episode_id"] = patient_request.get("episode_id")
    response["patient_id"] = patient_request.get("patient_id")

    return response
