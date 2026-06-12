"""
Micro-simulation engine for within-host bacterial evolution.

Runs 12 discrete time steps per day, simulating:
- Replication with fitness-based selection
- Mutation (Gaussian noise, stress-induced)
- Horizontal Gene Transfer (HGT) between strains
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Tuple

import numpy as np
from funkybob import UniqueRandomNameGenerator
from funkybob import data as funkybob_data

from .config import SimulationConfig
from .models import FounderStrain
from .genome import (
    ABX_PROFILES,
    DOSE_MULTIPLIERS,
    NUM_GENES,
    GeneIndex,
    classify_genotype,
    compute_abx_survival,
    compute_fitness,
    compute_resistant_fraction,
    compute_severity,
    compute_transmissibility,
    create_resistant_genome,
    create_wild_type_genome,
)


def _next_unique_strain_name(
    existing_names: set[str],
    rng: np.random.Generator,
    preferred_name: str | None = None,
) -> str:
    """Generate a unique two-word strain name.

    Args:
        existing_names: Names already assigned in the current population.
        rng: Random number generator used to seed the name generator.
        preferred_name: Optional name to use when it is not already taken.

    Returns:
        Unique strain name and records it in ``existing_names``.
    """
    if preferred_name and preferred_name not in existing_names:
        existing_names.add(preferred_name)
        return preferred_name

    while True:
        seed = int(rng.integers(0, 2**31 - 1))
        generator = UniqueRandomNameGenerator(
            members=2,
            seed=seed,
            names=funkybob_data.ADJECTIVES,
            adjectives=funkybob_data.ADJECTIVES,
        )
        name = generator[0]
        first, second = name.split("_", maxsplit=1)
        if first != second and name not in existing_names:
            existing_names.add(name)
            return name


def _extend_strain_names(
    names: list[str],
    count: int,
    rng: np.random.Generator,
    preferred_names: list[str] | None = None,
) -> list[str]:
    """Append unique strain names to an existing name list.

    Args:
        names: Existing strain names.
        count: Number of names to append.
        rng: Random number generator used to create names.
        preferred_names: Optional preferred names for the appended positions.

    Returns:
        New list containing the original and appended names.
    """
    updated = list(names)
    existing = set(updated)
    preferred_names = preferred_names or []
    for idx in range(count):
        preferred = preferred_names[idx] if idx < len(preferred_names) else None
        updated.append(_next_unique_strain_name(existing, rng, preferred_name=preferred))
    return updated


def _format_strain_id(namespace: str, serial: int) -> str:
    """Format a stable strain identifier.

    Args:
        namespace: Episode or population namespace.
        serial: Monotonic strain serial within the namespace.

    Returns:
        Stable strain identifier string.
    """
    return f"{namespace}:strain_{serial:06d}"


def _take_next_strain_id(namespace: str, next_serial: int) -> tuple[str, int]:
    """Allocate the next strain identifier and increment the serial.

    Args:
        namespace: Episode or population namespace.
        next_serial: Next available strain serial.

    Returns:
        Tuple of the allocated identifier and the following serial.
    """
    return _format_strain_id(namespace, next_serial), next_serial + 1


def _seed_founder_id(namespace: str, slot_index: int) -> str:
    """Format a fallback founder identifier for generated initial strains.

    Args:
        namespace: Episode or population namespace.
        slot_index: Initial strain slot index.

    Returns:
        Stable founder identifier string.
    """
    return f"{namespace}:founder_{slot_index:04d}"


def _generate_founder_genome(
    genotype: str,
    rng: np.random.Generator,
    noise_std: float,
) -> np.ndarray:
    """Generate a founder genome that remains in the requested genotype class.

    Args:
        genotype: Target genotype class, such as ``S``, ``R1``, ``R2``, or ``R3``.
        rng: Random number generator used for genome noise.
        noise_std: Standard deviation of founder gene noise.

    Returns:
        Genome vector classified as the requested genotype when possible.
    """
    base = (
        create_wild_type_genome() if genotype == "S" else _create_seed_genome_for_genotype(genotype)
    )
    if noise_std <= 0.0:
        return base.astype(np.float32, copy=True)

    # Resample (up to 16x) until the noised genome still classifies as the target
    # genotype, so founders stay within their intended resistance class.
    for _ in range(16):
        candidate = base.copy()
        candidate += rng.normal(0.0, noise_std, NUM_GENES).astype(np.float32)
        candidate = np.clip(candidate, 0.0, 1.0)
        if classify_genotype(candidate) == genotype:
            return candidate

    return base.astype(np.float32, copy=True)


def build_founder_pool(config: SimulationConfig) -> list[FounderStrain]:
    """Create a deterministic founder library from the active micro config.

    Args:
        config: Complete micro simulation configuration.

    Returns:
        List of founder strains. Returns an empty list when founder pooling is disabled.
    """
    size = max(0, int(config.founder_pool_size))
    if size == 0:
        return []

    forced_genotypes = ["S", "R1", "R2", "R3"]
    random_genotypes = np.array(["S", "R1", "R2", "R3"], dtype=object)
    random_weights = np.array([0.55, 0.20, 0.18, 0.07], dtype=np.float64)

    founders: list[FounderStrain] = []
    existing_names: set[str] = set()
    for offset in range(size):
        seed = int(config.founder_pool_seed) + offset
        rng = np.random.default_rng(seed)
        if offset < len(forced_genotypes):
            target_genotype = forced_genotypes[offset]
        else:
            target_genotype = str(rng.choice(random_genotypes, p=random_weights))

        genome = _generate_founder_genome(
            target_genotype,
            rng,
            noise_std=float(config.founder_pool_gene_noise_std),
        )
        founders.append(
            FounderStrain(
                founder_id=f"founder_{offset + 1:04d}",
                founder_name=_next_unique_strain_name(existing_names, rng),
                genotype=classify_genotype(genome),
                genome=genome.astype(np.float32, copy=True),
            )
        )

    return founders


def _select_founders_for_genotype(
    founders: list[FounderStrain],
    genotype: str,
    count: int,
    rng: np.random.Generator,
    preferred_name: str | None = None,
) -> list[FounderStrain]:
    """Select founder strains matching a genotype for initial population seeding.

    Args:
        founders: Candidate founder library.
        genotype: Required genotype class.
        count: Number of founders to select.
        rng: Random number generator used for sampling.
        preferred_name: Optional founder name to prioritize.

    Returns:
        Selected founders, sampling with replacement only when the eligible pool is small.
    """
    if count <= 0:
        return []

    eligible = [founder for founder in founders if founder.genotype == genotype]
    if not eligible:
        return []

    selected: list[FounderStrain] = []
    if preferred_name:
        preferred = next(
            (founder for founder in eligible if founder.founder_name == preferred_name),
            None,
        )
        if preferred is not None:
            selected.append(preferred)
            eligible = [
                founder for founder in eligible if founder.founder_id != preferred.founder_id
            ]

    remaining = count - len(selected)
    if remaining <= 0:
        return selected

    if len(eligible) >= remaining:
        indices = np.asarray(rng.choice(len(eligible), size=remaining, replace=False)).reshape(-1)
        selected.extend(eligible[int(idx)] for idx in indices.tolist())
        return selected

    selected.extend(eligible)
    remaining = count - len(selected)
    if remaining <= 0:
        return selected

    refill_pool = [founder for founder in founders if founder.genotype == genotype]
    indices = np.asarray(rng.choice(len(refill_pool), size=remaining, replace=True)).reshape(-1)
    selected.extend(refill_pool[int(idx)] for idx in indices.tolist())
    return selected


@dataclass
class StrainPopulation:
    """Discrete within-host bacterial population represented as strain arrays.

    Attributes:
        genomes: Genome matrix with shape ``(n_strains, NUM_GENES)``.
        populations: Population size for each strain.
        lineage_ages: Cumulative lineage age per strain.
        damage_loads: Accumulated stress or damage load per strain.
        strain_names: Human-readable strain names.
        strain_ids: Stable strain identifiers.
        parent_ids: Parent strain identifiers for mutant or recombinant strains.
        donor_ids: Donor strain identifiers for HGT-derived strains.
        founder_ids: Founder identifiers inherited by each strain lineage.
        next_strain_serial: Next serial used when allocating strain identifiers.
        strain_namespace: Namespace prefix used for generated strain identifiers.
    """

    genomes: np.ndarray  # Shape (n_strains, NUM_GENES)
    populations: np.ndarray  # Shape (n_strains,) - population sizes
    lineage_ages: np.ndarray | None = None  # Shape (n_strains,) - cumulative lineage age
    damage_loads: np.ndarray | None = None  # Shape (n_strains,) - accumulated stress / damage
    strain_names: list[str] | None = None
    strain_ids: list[str] | None = None
    parent_ids: list[str] | None = None
    donor_ids: list[str] | None = None
    founder_ids: list[str] | None = None
    next_strain_serial: int = 0
    strain_namespace: str = "population"

    def __post_init__(self):
        """Normalize optional arrays and validate strain-aligned fields."""
        assert self.genomes.ndim == 2
        assert self.genomes.shape[1] == NUM_GENES
        assert len(self.populations) == len(self.genomes)
        self.strain_namespace = str(self.strain_namespace or "population")
        if self.lineage_ages is None:
            self.lineage_ages = np.zeros(len(self.populations), dtype=np.float64)
        if self.damage_loads is None:
            self.damage_loads = np.zeros(len(self.populations), dtype=np.float64)
        self.lineage_ages = np.asarray(self.lineage_ages, dtype=np.float64)
        self.damage_loads = np.asarray(self.damage_loads, dtype=np.float64)
        assert len(self.lineage_ages) == len(self.genomes)
        assert len(self.damage_loads) == len(self.genomes)
        if self.strain_names is None:
            self.strain_names = [f"strain_{idx:04d}" for idx in range(len(self.populations))]
        else:
            self.strain_names = [str(name) for name in self.strain_names]
        assert len(self.strain_names) == len(self.genomes)
        if self.strain_ids is None:
            self.strain_ids = [
                _format_strain_id(self.strain_namespace, idx)
                for idx in range(len(self.populations))
            ]
        else:
            self.strain_ids = [str(value) for value in self.strain_ids]
        if self.parent_ids is None:
            self.parent_ids = ["" for _ in range(len(self.populations))]
        else:
            self.parent_ids = [str(value) for value in self.parent_ids]
        if self.donor_ids is None:
            self.donor_ids = ["" for _ in range(len(self.populations))]
        else:
            self.donor_ids = [str(value) for value in self.donor_ids]
        if self.founder_ids is None:
            self.founder_ids = self.strain_ids.copy()
        else:
            self.founder_ids = [str(value) for value in self.founder_ids]
        assert len(self.strain_ids) == len(self.genomes)
        assert len(self.parent_ids) == len(self.genomes)
        assert len(self.donor_ids) == len(self.genomes)
        assert len(self.founder_ids) == len(self.genomes)
        self.next_strain_serial = max(int(self.next_strain_serial), len(self.strain_ids))

    @property
    def n_strains(self) -> int:
        """Return the number of active strains.

        Returns:
            Number of population entries.
        """
        return len(self.populations)

    @property
    def total_population(self) -> float:
        """Return the total within-host bacterial population.

        Returns:
            Sum of all strain populations.
        """
        return float(np.sum(self.populations))

    def clone(self) -> StrainPopulation:
        """Create a deep copy of the strain population.

        Returns:
            Independent population object with copied arrays and metadata lists.
        """
        return StrainPopulation(
            genomes=self.genomes.copy(),
            populations=self.populations.copy(),
            lineage_ages=self.lineage_ages.copy(),
            damage_loads=self.damage_loads.copy(),
            strain_names=self.strain_names.copy(),
            strain_ids=self.strain_ids.copy(),
            parent_ids=self.parent_ids.copy(),
            donor_ids=self.donor_ids.copy(),
            founder_ids=self.founder_ids.copy(),
            next_strain_serial=self.next_strain_serial,
            strain_namespace=self.strain_namespace,
        )

    @classmethod
    def create_initial(
        cls,
        resistant_fraction: float = 0.0,
        dominant_genotype: str = "S",
        initial_population: float = 1e6,
        n_susceptible_strains: int = 3,
        n_resistant_strains: int = 2,
        rng: np.random.Generator = None,
        dominant_strain_name: str | None = None,
        founder_pool: list[FounderStrain] | None = None,
        strain_namespace: str = "population",
        seed_genome: np.ndarray | None = None,
    ) -> StrainPopulation:
        """Create an initial strain population for a new episode.

        When ``seed_genome`` is provided, all strains are initialized as variations
        around that transmitted genome. Otherwise, strains are drawn from the founder
        pool when available or generated from susceptible/resistant templates.

        Args:
            resistant_fraction: Initial fraction assigned to resistant strains.
            dominant_genotype: Target dominant genotype class.
            initial_population: Total initial bacterial population size.
            n_susceptible_strains: Number of susceptible initial strains.
            n_resistant_strains: Number of resistant initial strains.
            rng: Random number generator. A default generator is created when omitted.
            dominant_strain_name: Optional preferred name for the dominant initial strain.
            founder_pool: Optional reusable founder strain library.
            strain_namespace: Namespace used for generated strain identifiers.
            seed_genome: Optional transmitted genome used to seed the population.

        Returns:
            Initialized strain population with strain IDs and lineage metadata.
        """
        if rng is None:
            rng = np.random.default_rng()

        strains = []
        pops = []
        names: list[str] = []
        strain_ids: list[str] = []
        parent_ids: list[str] = []
        donor_ids: list[str] = []
        founder_ids: list[str] = []
        existing_names: set[str] = set()
        next_serial = 0
        founder_pool = founder_pool or []

        # Susceptible strains
        sus_pop = initial_population * (1 - resistant_fraction)

        if seed_genome is None:
            susceptible_founders = _select_founders_for_genotype(
                founder_pool,
                "S",
                n_susceptible_strains,
                rng,
                preferred_name=dominant_strain_name if dominant_genotype == "S" else None,
            )

        for i in range(n_susceptible_strains):
            if seed_genome is not None:
                genome = (
                    seed_genome.copy().astype(np.float32)
                    if i == 0
                    else np.clip(
                        seed_genome + rng.normal(0, 0.02, NUM_GENES).astype(np.float32),
                        0.0,
                        1.0,
                    )
                )
                preferred_name = (
                    dominant_strain_name if dominant_genotype == "S" and i == 0 else None
                )
                founder_id = _seed_founder_id(strain_namespace, len(founder_ids))
            else:
                founder = susceptible_founders[i] if i < len(susceptible_founders) else None
                if founder is not None:
                    genome = founder.genome.copy()
                    preferred_name = founder.founder_name
                    founder_id = founder.founder_id
                else:
                    genome = create_wild_type_genome()
                    genome += rng.normal(0, 0.02, NUM_GENES).astype(np.float32)
                    genome = np.clip(genome, 0.0, 1.0)
                    preferred_name = (
                        dominant_strain_name if dominant_genotype == "S" and i == 0 else None
                    )
                    founder_id = _seed_founder_id(strain_namespace, len(founder_ids))
            strains.append(genome)
            pops.append(sus_pop / n_susceptible_strains)
            names.append(
                _next_unique_strain_name(existing_names, rng, preferred_name=preferred_name)
            )
            strain_id, next_serial = _take_next_strain_id(strain_namespace, next_serial)
            strain_ids.append(strain_id)
            parent_ids.append("")
            donor_ids.append("")
            founder_ids.append(founder_id)

        # Resistant strains (if any)
        if resistant_fraction > 0 and n_resistant_strains > 0:
            res_pop = initial_population * resistant_fraction
            resistant_target = dominant_genotype if dominant_genotype != "S" else "R1"

            if seed_genome is None:
                resistant_founders = _select_founders_for_genotype(
                    founder_pool,
                    resistant_target,
                    n_resistant_strains,
                    rng,
                    preferred_name=dominant_strain_name if dominant_genotype != "S" else None,
                )

            for i in range(n_resistant_strains):
                if seed_genome is not None:
                    genome = np.clip(
                        seed_genome + rng.normal(0, 0.02, NUM_GENES).astype(np.float32),
                        0.0,
                        1.0,
                    )
                    preferred_name = (
                        dominant_strain_name if i == 0 and dominant_genotype != "S" else None
                    )
                    founder_id = _seed_founder_id(strain_namespace, len(founder_ids))
                else:
                    founder = resistant_founders[i] if i < len(resistant_founders) else None
                    if founder is not None:
                        genome = founder.genome.copy()
                        preferred_name = founder.founder_name
                        founder_id = founder.founder_id
                    else:
                        genome = _create_seed_genome_for_genotype(dominant_genotype)
                        genome += rng.normal(0, 0.02, NUM_GENES).astype(np.float32)
                        genome = np.clip(genome, 0.0, 1.0)
                        preferred_name = (
                            dominant_strain_name if i == 0 and dominant_genotype != "S" else None
                        )
                        founder_id = _seed_founder_id(strain_namespace, len(founder_ids))
                strains.append(genome)
                pops.append(res_pop / n_resistant_strains)
                names.append(
                    _next_unique_strain_name(existing_names, rng, preferred_name=preferred_name)
                )
                strain_id, next_serial = _take_next_strain_id(strain_namespace, next_serial)
                strain_ids.append(strain_id)
                parent_ids.append("")
                donor_ids.append("")
                founder_ids.append(founder_id)

        genomes = np.array(strains, dtype=np.float32)
        populations = np.array(pops, dtype=np.float64)
        lineage_ages = np.zeros(len(populations), dtype=np.float64)
        damage_loads = np.zeros(len(populations), dtype=np.float64)

        return cls(
            genomes=genomes,
            populations=populations,
            lineage_ages=lineage_ages,
            damage_loads=damage_loads,
            strain_names=names,
            strain_ids=strain_ids,
            parent_ids=parent_ids,
            donor_ids=donor_ids,
            founder_ids=founder_ids,
            next_strain_serial=next_serial,
            strain_namespace=strain_namespace,
        )


def _create_seed_genome_for_genotype(dominant_genotype: str) -> np.ndarray:
    """Create a seed genome that classifies as the requested genotype.

    Args:
        dominant_genotype: Target genotype class.

    Returns:
        Genome vector for the requested genotype, falling back to low resistance for
        unknown resistant labels.
    """
    genome = create_wild_type_genome()

    if dominant_genotype == "R3":
        genome[GeneIndex.EFFLUX_PUMPS] = 0.95
        genome[GeneIndex.TARGET_MODIFICATION] = 0.85
        genome[GeneIndex.PERMEABILITY_REDUCTION] = 0.75
        genome[GeneIndex.METABOLIC_OPTIMIZATION] = 0.65
        genome[GeneIndex.STRESS_RESPONSE] = 0.65
        genome[GeneIndex.DAMAGE_TOLERANCE] = 0.55
        genome[GeneIndex.DORMANCY_PROPENSITY] = 0.35
        genome[GeneIndex.GROWTH_BASE] = 0.62
        return genome

    if dominant_genotype == "R2":
        genome[GeneIndex.EFFLUX_PUMPS] = 0.60
        genome[GeneIndex.TARGET_MODIFICATION] = 0.50
        genome[GeneIndex.PERMEABILITY_REDUCTION] = 0.40
        genome[GeneIndex.METABOLIC_OPTIMIZATION] = 0.45
        genome[GeneIndex.STRESS_RESPONSE] = 0.45
        genome[GeneIndex.DAMAGE_TOLERANCE] = 0.35
        genome[GeneIndex.DORMANCY_PROPENSITY] = 0.22
        genome[GeneIndex.GROWTH_BASE] = 0.68
        return genome

    resistance_level = 0.35 if dominant_genotype == "R1" else 0.30
    return create_resistant_genome(resistance_level)


def mutate_population(
    population: StrainPopulation,
    config: SimulationConfig,
    abx_stress: float,
    rng: np.random.Generator,
) -> StrainPopulation:
    """Apply mutation events and create mutant strains.

    Args:
        population: Current strain population.
        config: Complete micro simulation configuration.
        abx_stress: Antibiotic stress level in ``[0, 1]`` used to boost mutation rate.
        rng: Random number generator used for event sampling and mutation effects.

    Returns:
        Updated population, possibly containing new mutant strains and adjusted parent
        populations.
    """
    genomes = population.genomes.copy()
    populations = population.populations.copy()
    lineage_ages = population.lineage_ages.copy()
    damage_loads = population.damage_loads.copy()
    strain_names = population.strain_names.copy()
    strain_ids = population.strain_ids.copy()
    parent_ids = population.parent_ids.copy()
    donor_ids = population.donor_ids.copy()
    founder_ids = population.founder_ids.copy()
    next_serial = population.next_strain_serial

    # Stress-induced mutation rate increase
    stress_mult = 1.0 + abx_stress * (config.stress_mutation_boost - 1.0)
    effective_rate = config.base_mutation_rate * stress_mult

    # Mutation rate modifier from genome
    mutation_mods = genomes[:, GeneIndex.MUTATION_RATE_MODIFIER]
    strain_rates = effective_rate * (0.5 + mutation_mods)

    new_strains = []
    new_pops = []
    new_ages = []
    new_damage = []
    new_ids: list[str] = []
    new_parent_ids: list[str] = []
    new_donor_ids: list[str] = []
    new_founder_ids: list[str] = []

    for i in range(len(genomes)):
        if populations[i] < 1:
            continue

        # Number of mutations this step (Poisson)
        n_mutations = rng.poisson(strain_rates[i] * NUM_GENES)

        if n_mutations > 0 and populations[i] > 100:
            # Create mutant strain
            mutant = genomes[i].copy()

            # Select genes to mutate
            genes_to_mutate = rng.choice(NUM_GENES, size=min(n_mutations, NUM_GENES), replace=False)

            for gene_idx in genes_to_mutate:
                delta = rng.normal(0, config.mutation_std)
                mutant[gene_idx] = np.clip(mutant[gene_idx] + delta, 0.0, 1.0)

            # Transfer a configurable founder fraction so mutants are not born below pruning.
            transfer_fraction = min(0.10, config.mutant_transfer_fraction * n_mutations)
            transfer_pop = populations[i] * transfer_fraction
            populations[i] -= transfer_pop

            new_strains.append(mutant)
            new_pops.append(transfer_pop)
            new_ages.append(max(0.0, lineage_ages[i] * 0.5))
            new_damage.append(max(0.0, damage_loads[i] * 0.7))
            strain_id, next_serial = _take_next_strain_id(population.strain_namespace, next_serial)
            new_ids.append(strain_id)
            new_parent_ids.append(strain_ids[i])
            new_donor_ids.append("")
            new_founder_ids.append(founder_ids[i])

    # Add new strains
    if new_strains:
        genomes = np.vstack([genomes, np.array(new_strains, dtype=np.float32)])
        populations = np.concatenate([populations, np.array(new_pops)])
        lineage_ages = np.concatenate([lineage_ages, np.array(new_ages, dtype=np.float64)])
        damage_loads = np.concatenate([damage_loads, np.array(new_damage, dtype=np.float64)])
        strain_names = _extend_strain_names(strain_names, len(new_strains), rng)
        strain_ids.extend(new_ids)
        parent_ids.extend(new_parent_ids)
        donor_ids.extend(new_donor_ids)
        founder_ids.extend(new_founder_ids)

    return StrainPopulation(
        genomes=genomes,
        populations=populations,
        lineage_ages=lineage_ages,
        damage_loads=damage_loads,
        strain_names=strain_names,
        strain_ids=strain_ids,
        parent_ids=parent_ids,
        donor_ids=donor_ids,
        founder_ids=founder_ids,
        next_strain_serial=next_serial,
        strain_namespace=population.strain_namespace,
    )


def horizontal_gene_transfer(
    population: StrainPopulation, config: SimulationConfig, rng: np.random.Generator
) -> StrainPopulation:
    """Simulate horizontal gene transfer between strains.

    Args:
        population: Current strain population.
        config: Complete micro simulation configuration.
        rng: Random number generator used for donor selection and gene blending.

    Returns:
        Updated population, possibly containing recombinant strains.
    """
    if population.n_strains < 2:
        return population

    genomes = population.genomes.copy()
    populations = population.populations.copy()
    lineage_ages = population.lineage_ages.copy()
    damage_loads = population.damage_loads.copy()
    strain_names = population.strain_names.copy()
    strain_ids = population.strain_ids.copy()
    parent_ids = population.parent_ids.copy()
    donor_ids = population.donor_ids.copy()
    founder_ids = population.founder_ids.copy()
    next_serial = population.next_strain_serial

    # HGT probability based on competence
    competence = genomes[:, GeneIndex.HGT_COMPETENCE]

    # Genes that can be horizontally transferred (mainly resistance)
    transferable_genes = [
        GeneIndex.EFFLUX_PUMPS,
        GeneIndex.TARGET_MODIFICATION,
        GeneIndex.PERMEABILITY_REDUCTION,
        GeneIndex.METABOLIC_OPTIMIZATION,
        GeneIndex.DORMANCY_PROPENSITY,
        GeneIndex.STRESS_RESPONSE,
        GeneIndex.DAMAGE_TOLERANCE,
    ]

    new_strains = []
    new_pops = []
    new_ages = []
    new_damage = []
    new_ids: list[str] = []
    new_parent_ids: list[str] = []
    new_donor_ids: list[str] = []
    new_founder_ids: list[str] = []

    for i in range(population.n_strains):
        if populations[i] < 1000:
            continue

        hgt_prob = config.base_hgt_rate * (0.5 + competence[i])

        if rng.random() < hgt_prob:
            # Select donor strain (weighted by population)
            donor_weights = populations.copy()
            donor_weights[i] = 0  # Can't be own donor
            if np.sum(donor_weights) == 0:
                continue
            donor_weights /= np.sum(donor_weights)

            donor_idx = rng.choice(population.n_strains, p=donor_weights)

            # Create recombinant
            recombinant = genomes[i].copy()

            for gene_idx in transferable_genes:
                if rng.random() < config.hgt_gene_transfer_prob:
                    # Blend genes (partial transfer)
                    blend = 0.3 + rng.random() * 0.4
                    recombinant[gene_idx] = (
                        recombinant[gene_idx] * (1 - blend) + genomes[donor_idx, gene_idx] * blend
                    )

            # Small fraction becomes recombinant
            transfer_pop = populations[i] * 0.005
            populations[i] -= transfer_pop

            new_strains.append(recombinant)
            new_pops.append(transfer_pop)
            new_ages.append(max(0.0, min(lineage_ages[i], lineage_ages[donor_idx]) * 0.7))
            new_damage.append(max(0.0, 0.5 * (damage_loads[i] + damage_loads[donor_idx]) * 0.8))
            strain_id, next_serial = _take_next_strain_id(population.strain_namespace, next_serial)
            new_ids.append(strain_id)
            new_parent_ids.append(strain_ids[i])
            new_donor_ids.append(strain_ids[donor_idx])
            new_founder_ids.append(founder_ids[i])

    if new_strains:
        genomes = np.vstack([genomes, np.array(new_strains, dtype=np.float32)])
        populations = np.concatenate([populations, np.array(new_pops)])
        lineage_ages = np.concatenate([lineage_ages, np.array(new_ages, dtype=np.float64)])
        damage_loads = np.concatenate([damage_loads, np.array(new_damage, dtype=np.float64)])
        strain_names = _extend_strain_names(strain_names, len(new_strains), rng)
        strain_ids.extend(new_ids)
        parent_ids.extend(new_parent_ids)
        donor_ids.extend(new_donor_ids)
        founder_ids.extend(new_founder_ids)

    return StrainPopulation(
        genomes=genomes,
        populations=populations,
        lineage_ages=lineage_ages,
        damage_loads=damage_loads,
        strain_names=strain_names,
        strain_ids=strain_ids,
        parent_ids=parent_ids,
        donor_ids=donor_ids,
        founder_ids=founder_ids,
        next_strain_serial=next_serial,
        strain_namespace=population.strain_namespace,
    )


def apply_demographic_stochasticity(
    populations: np.ndarray,
    config: SimulationConfig,
    rng: np.random.Generator,
) -> np.ndarray:
    """Apply finite-population demographic noise to strain sizes.

    Args:
        populations: Expected strain populations after deterministic update.
        config: Complete micro simulation configuration.
        rng: Random number generator used for Poisson or normal sampling.

    Returns:
        Realized non-negative strain populations.
    """
    realized = np.zeros_like(populations, dtype=np.float64)

    for i, expected in enumerate(populations):
        if expected <= 0.0:
            continue
        if expected <= config.stochastic_threshold:
            realized[i] = float(rng.poisson(expected))
            continue

        sd = np.sqrt(expected) * config.stochastic_noise_scale
        realized[i] = float(max(0.0, rng.normal(expected, sd)))

    return realized


def selection_step(
    population: StrainPopulation,
    config: SimulationConfig,
    abx_class: str,
    dose_level: str,
    adherence: float,
    immune_strength: float,
    rng: np.random.Generator,
) -> StrainPopulation:
    """Advance one selection/growth/death step.

    Args:
        population: Current strain population.
        config: Complete micro simulation configuration.
        abx_class: Active antibiotic class, or ``none``.
        dose_level: Active antibiotic dose label.
        adherence: Effective adherence multiplier.
        immune_strength: Host immune-strength multiplier.
        rng: Random number generator used for demographic stochasticity.

    Returns:
        Population after density-regulated growth, explicit ABX kill, lifecycle turnover,
        and demographic stochasticity.
    """
    if population.n_strains == 0:
        return population

    genomes = population.genomes
    populations = population.populations.copy()
    lineage_ages = population.lineage_ages.copy()
    damage_loads = population.damage_loads.copy()

    # Baseline host fitness is separate from antibiotic killing. Antibiotics
    # create directional selection and an explicit kill term below.
    baseline_fitness = compute_fitness(
        genomes,
        abx_class="none",
        dose_level=dose_level,
        adherence=adherence,
        immune_strength=immune_strength,
    )
    abx_survival = compute_abx_survival(
        genomes,
        abx_class=abx_class,
        dose_level=dose_level,
        adherence=adherence,
    )
    selection_fitness = np.clip(baseline_fitness * abx_survival, 0.001, 1.0)

    # Relative fitness affects growth rate. Under ABX this favors resistant tails;
    # without ABX it reduces to baseline host fitness.
    mean_fitness = np.mean(selection_fitness) if len(selection_fitness) > 0 else 0.5
    relative_fitness = selection_fitness / (mean_fitness + 1e-6)

    # Antibiotic pressure turns small fitness differences into life-or-death selection.
    profile = ABX_PROFILES.get(abx_class, ABX_PROFILES["none"])
    dose_multiplier = DOSE_MULTIPLIERS.get(dose_level, 1.0)
    abx_pressure = float(np.clip(profile.base_kill_rate * dose_multiplier * adherence, 0.0, 1.0))
    effective_selection_strength = config.selection_strength * (
        1.0 + abx_pressure * config.abx_selection_pressure_multiplier
    )
    selection_factor = np.power(relative_fitness, effective_selection_strength)

    # Population dynamics
    growth = config.growth_rate_per_step * selection_factor * selection_fitness
    metabolic_support = np.clip(genomes[:, GeneIndex.METABOLIC_OPTIMIZATION], 0.0, 1.0)
    repair_gene = np.clip(genomes[:, GeneIndex.DNA_REPAIR], 0.0, 1.0)
    dormancy_gene = np.clip(genomes[:, GeneIndex.DORMANCY_PROPENSITY], 0.0, 1.0)
    stress_gene = np.clip(genomes[:, GeneIndex.STRESS_RESPONSE], 0.0, 1.0)
    tolerance_gene = np.clip(genomes[:, GeneIndex.DAMAGE_TOLERANCE], 0.0, 1.0)

    repair_dormancy_combo = repair_gene * dormancy_gene
    stress_tolerance_combo = stress_gene * tolerance_gene

    repair_capacity = np.clip(
        0.65 * repair_gene
        + 0.20 * metabolic_support
        + 0.15 * stress_gene
        + config.synergy_repair_dormancy_bonus * repair_dormancy_combo,
        0.0,
        1.5,
    )
    dormancy_capacity = np.clip(
        0.75 * dormancy_gene + 0.15 * stress_gene + 0.10 * genomes[:, GeneIndex.STEALTH],
        0.0,
        1.2,
    )
    tolerance_capacity = np.clip(
        0.70 * tolerance_gene
        + 0.20 * stress_gene
        + 0.10 * metabolic_support
        + config.synergy_stress_tolerance_bonus * stress_tolerance_combo,
        0.0,
        1.5,
    )

    environmental_stress = 1.0 - selection_fitness
    active_dormancy = dormancy_capacity * (0.2 + 0.8 * environmental_stress)

    growth *= np.clip(1.0 - config.dormancy_growth_penalty * active_dormancy, 0.2, 1.0)

    # Logistic density regulation: carrying_capacity slows growth as the
    # within-host population fills the niche and reverses growth above capacity.
    total_population = float(np.sum(populations))
    if config.carrying_capacity > 0.0:
        density_factor = 1.0 - (total_population / config.carrying_capacity)
        growth *= np.clip(density_factor, -5.0, 1.0)

    replication_pressure = np.maximum(growth, 0.0) * (1.0 - 0.7 * active_dormancy)
    damage_gain = (
        config.base_damage_per_step
        + config.replication_damage_factor * replication_pressure
        + config.stress_damage_factor * environmental_stress
    )
    damage_gain *= np.clip(1.0 - 0.45 * repair_capacity - 0.20 * tolerance_capacity, 0.10, 1.0)

    repair = config.repair_rate_per_step * (
        0.15 + repair_capacity + 0.35 * active_dormancy + 0.25 * stress_tolerance_combo
    )
    damage_loads = np.clip(damage_loads + damage_gain - repair, 0.0, config.max_damage_load)

    lineage_ages += 1.0 + replication_pressure
    age_pressure = lineage_ages / config.lifecycle_half_life_steps
    damage_pressure = damage_loads / config.max_damage_load
    turnover_pressure = (
        config.age_mortality_scale * age_pressure + config.damage_mortality_scale * damage_pressure
    )
    turnover_pressure *= np.clip(
        1.0
        - 0.35 * repair_capacity
        - 0.20 * active_dormancy
        - 0.25 * tolerance_capacity
        - 0.10 * repair_dormancy_combo,
        0.08,
        1.0,
    )

    baseline_death = config.death_rate_per_step * (
        1.0 / (baseline_fitness + config.death_fitness_floor)
    )
    antibiotic_kill = config.antibiotic_kill_scale * (1.0 - abx_survival)
    death = baseline_death + antibiotic_kill + turnover_pressure

    # Net growth uses exponential compounding so rates behave like continuous hazards.
    net_growth = growth - death
    populations = populations * np.exp(np.clip(net_growth, -50.0, 50.0))
    populations = apply_demographic_stochasticity(populations, config, rng)

    return StrainPopulation(
        genomes=genomes,
        populations=populations,
        lineage_ages=lineage_ages,
        damage_loads=damage_loads,
        strain_names=population.strain_names.copy(),
        strain_ids=population.strain_ids.copy(),
        parent_ids=population.parent_ids.copy(),
        donor_ids=population.donor_ids.copy(),
        founder_ids=population.founder_ids.copy(),
        next_strain_serial=population.next_strain_serial,
        strain_namespace=population.strain_namespace,
    )


def consolidate_strains(population: StrainPopulation, config: SimulationConfig) -> StrainPopulation:
    """Prune small strains and enforce the maximum tracked strain count.

    Args:
        population: Current strain population.
        config: Complete micro simulation configuration.

    Returns:
        Population containing only strains above the prune threshold and, when needed,
        the largest ``max_strains`` strains.
    """
    # Remove strains below threshold
    mask = population.populations > config.strain_prune_threshold

    if not np.any(mask):
        return StrainPopulation(
            genomes=np.empty((0, NUM_GENES), dtype=np.float32),
            populations=np.empty(0, dtype=np.float64),
            lineage_ages=np.empty(0, dtype=np.float64),
            damage_loads=np.empty(0, dtype=np.float64),
            strain_names=[],
            strain_ids=[],
            parent_ids=[],
            donor_ids=[],
            founder_ids=[],
            next_strain_serial=population.next_strain_serial,
            strain_namespace=population.strain_namespace,
        )

    genomes = population.genomes[mask]
    populations = population.populations[mask]
    lineage_ages = population.lineage_ages[mask]
    damage_loads = population.damage_loads[mask]
    strain_names = [name for name, keep in zip(population.strain_names, mask, strict=False) if keep]
    strain_ids = [value for value, keep in zip(population.strain_ids, mask, strict=False) if keep]
    parent_ids = [value for value, keep in zip(population.parent_ids, mask, strict=False) if keep]
    donor_ids = [value for value, keep in zip(population.donor_ids, mask, strict=False) if keep]
    founder_ids = [value for value, keep in zip(population.founder_ids, mask, strict=False) if keep]

    # If too many strains, keep the largest ones
    if len(populations) > config.max_strains:
        indices = np.argsort(populations)[-config.max_strains :]
        genomes = genomes[indices]
        populations = populations[indices]
        lineage_ages = lineage_ages[indices]
        damage_loads = damage_loads[indices]
        strain_names = [strain_names[idx] for idx in indices]
        strain_ids = [strain_ids[idx] for idx in indices]
        parent_ids = [parent_ids[idx] for idx in indices]
        donor_ids = [donor_ids[idx] for idx in indices]
        founder_ids = [founder_ids[idx] for idx in indices]

    return StrainPopulation(
        genomes=genomes,
        populations=populations,
        lineage_ages=lineage_ages,
        damage_loads=damage_loads,
        strain_names=strain_names,
        strain_ids=strain_ids,
        parent_ids=parent_ids,
        donor_ids=donor_ids,
        founder_ids=founder_ids,
        next_strain_serial=population.next_strain_serial,
        strain_namespace=population.strain_namespace,
    )


def simulate_day(
    population: StrainPopulation,
    abx_class: str,
    dose_level: str,
    adherence: float,
    immune_strength: float,
    config: SimulationConfig,
    seed: int = None,
) -> Tuple[StrainPopulation, Dict[str, Any]]:
    """Run one macro day of within-host micro dynamics.

    Args:
        population: Starting strain population.
        abx_class: Antibiotic class administered during the day, or ``none``.
        dose_level: Antibiotic dose label.
        adherence: Effective adherence multiplier.
        immune_strength: Host immune-strength multiplier.
        config: Complete micro simulation configuration.
        seed: Optional random seed for reproducible within-day dynamics.

    Returns:
        Tuple containing the final population and a step-history dictionary.
    """
    rng = np.random.default_rng(seed)

    # Compute ABX stress level
    profile = ABX_PROFILES.get(abx_class, ABX_PROFILES["none"])
    abx_stress = profile.base_kill_rate * adherence if profile.base_kill_rate > 0 else 0.0

    pop = population.clone()
    history = []

    for step in range(config.steps_per_day):
        # 1. Selection (growth/death based on fitness)
        pop = selection_step(pop, config, abx_class, dose_level, adherence, immune_strength, rng)

        # 2. Mutation
        pop = mutate_population(pop, config, abx_stress, rng)

        # 3. HGT (less frequent)
        if step % 3 == 0:  # Every 3rd step
            pop = horizontal_gene_transfer(pop, config, rng)

        # 4. Consolidate strains
        pop = consolidate_strains(pop, config)

        # Record state
        history.append(
            {
                "step": step,
                "total_pop": pop.total_population,
                "n_strains": pop.n_strains,
                "resistant_fraction": compute_resistant_fraction(pop.genomes, pop.populations),
            }
        )

    return pop, {"steps": history}


def compute_clearance_probability(
    population: StrainPopulation,
    immune_strength: float,
    config: SimulationConfig,
) -> float:
    """Compute daily bacterial clearance probability from micro state.

    Args:
        population: Current strain population.
        immune_strength: Host immune-strength multiplier.
        config: Complete micro simulation configuration.

    Returns:
        Probability that the carrier clears the episode on this macro day.
    """
    total_pop = population.total_population

    if total_pop <= 0:
        return 1.0

    if total_pop < config.min_population:
        return 1.0  # extinction floor: carriage cleared

    if total_pop < config.clearance_threshold:
        return 0.95  # very likely to clear

    # Logistic regime: 0.006 prefactor tuned so daily clearance matches the
    # macro-level p_clearance ~= 0.0039/day (~1/255 days mean MRSA carriage).
    ratio = total_pop / config.carrying_capacity
    base_prob = 0.006 * (1.0 - ratio)

    # Immune modulation
    immune_mult = immune_strength

    # Stealth of dominant strain reduces clearance
    if population.n_strains > 0:
        # Weight by population
        weights = population.populations / (np.sum(population.populations) + 1e-6)
        avg_stealth = np.sum(population.genomes[:, GeneIndex.STEALTH] * weights)
        stealth_effect = 1.0 - avg_stealth * 0.5
    else:
        stealth_effect = 1.0

    p_clear = base_prob * immune_mult * stealth_effect
    return float(np.clip(p_clear, 0.001, 0.95))


def get_dominant_strain(population: StrainPopulation) -> Tuple[np.ndarray, str, str]:
    """Get the dominant strain summary.

    Args:
        population: Current strain population.

    Returns:
        Tuple of dominant genome, genotype class, and strain name. Empty populations return
        a susceptible wild-type genome with an empty strain name.
    """
    if population.n_strains == 0:
        genome = create_wild_type_genome()
        return genome, "S", ""

    idx = np.argmax(population.populations)
    genome = population.genomes[idx]
    genotype = classify_genotype(genome)
    strain_name = population.strain_names[idx]

    return genome, genotype, strain_name


def population_to_response(
    population: StrainPopulation,
    immune_strength: float,
    config: SimulationConfig,
) -> Dict[str, Any]:
    """Convert a final micro population into a macro-compatible response.

    Args:
        population: Final strain population for the simulated day.
        immune_strength: Host immune-strength multiplier used for clearance probability.
        config: Complete micro simulation configuration.

    Returns:
        Response dictionary compatible with ``Patient.apply_micro_response``.
    """
    dominant_genome, dominant_genotype, dominant_strain_name = get_dominant_strain(population)
    resistant_fraction = compute_resistant_fraction(population.genomes, population.populations)

    # Compute derived effects from dominant strain
    transmissibility = float(compute_transmissibility(dominant_genome))
    severity = float(compute_severity(dominant_genome))
    p_clearance = compute_clearance_probability(population, immune_strength, config)

    return {
        "updated_state": {
            "resistant_fraction": resistant_fraction,
            "dominant_genotype": dominant_genotype,
            "dominant_strain_name": dominant_strain_name,
            "dominant_genome": dominant_genome.tolist(),
        },
        "derived_effects": {
            "relative_transmissibility": transmissibility,
            "severity_modifier": severity,
            "p_clearance": p_clearance,
        },
        "population_stats": {
            "total_population": population.total_population,
            "n_strains": population.n_strains,
            "dominant_strain_name": dominant_strain_name,
        },
    }
