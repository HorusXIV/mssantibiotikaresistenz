"""
Micro-simulation engine for within-host bacterial evolution.

Runs 12 discrete time steps per day, simulating:
- Replication with fitness-based selection
- Mutation (Gaussian noise, stress-induced)
- Horizontal Gene Transfer (HGT) between strains
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Dict, Any, Tuple

from .genome import (
    NUM_GENES,
    GeneIndex,
    create_wild_type_genome,
    create_resistant_genome,
    compute_fitness,
    compute_transmissibility,
    compute_lethality,
    compute_severity,
    classify_genotype,
    compute_resistant_fraction,
    ABX_PROFILES,
)


@dataclass
class SimulationConfig:
    """Configuration for micro-simulation."""

    steps_per_day: int = 12
    max_strains: int = 50  # Max distinct strains to track
    carrying_capacity: float = 1e9  # Max population size
    min_population: float = 1e3  # Below this, clearance likely
    clearance_threshold: float = 1e2  # Population for clearance

    # Mutation parameters
    base_mutation_rate: float = 0.01  # Base mutation per gene per step
    mutation_std: float = 0.05  # Gaussian std for mutations
    stress_mutation_boost: float = 3.0  # Mutation rate multiplier under ABX stress

    # HGT parameters
    base_hgt_rate: float = 0.02  # Base probability of HGT per step
    hgt_gene_transfer_prob: float = 0.3  # Prob of transferring each gene

    # Selection parameters
    selection_strength: float = 2.0  # Exponent for fitness-based selection

    # Population dynamics
    growth_rate_per_step: float = 0.3  # Max growth per step (before fitness)
    death_rate_per_step: float = 0.1  # Base death rate per step
    strain_prune_threshold: float = 1.0  # Drop numerically negligible strains

    # Lifecycle / turnover dynamics
    base_damage_per_step: float = 0.03  # Background wear per step
    replication_damage_factor: float = 0.25  # Extra wear from fast turnover
    stress_damage_factor: float = 0.20  # Extra wear under environmental stress
    repair_rate_per_step: float = 0.08  # Damage repair / maintenance capacity
    age_mortality_scale: float = 0.015  # How strongly lineage age raises death risk
    damage_mortality_scale: float = 0.20  # How strongly accumulated damage raises death risk
    lifecycle_half_life_steps: float = 36.0  # Age scale for senescence-like pressure
    max_damage_load: float = 5.0  # Saturation point for damage effects
    dormancy_growth_penalty: float = 0.25  # Persistence trades growth for survival
    synergy_repair_dormancy_bonus: float = 0.25  # Quiescent repair synergy
    synergy_stress_tolerance_bonus: float = 0.20  # Stress hardening synergy

    # Demographic stochasticity
    stochastic_threshold: float = 1e5  # Exact Poisson dynamics below this size
    stochastic_noise_scale: float = 1.0  # Scales demographic noise above threshold


@dataclass
class StrainPopulation:
    """
    Represents bacterial population as discrete strains with population counts.

    Uses vectorized NumPy arrays for efficient computation.
    """

    genomes: np.ndarray  # Shape (n_strains, NUM_GENES)
    populations: np.ndarray  # Shape (n_strains,) - population sizes
    lineage_ages: np.ndarray | None = None  # Shape (n_strains,) - cumulative lineage age
    damage_loads: np.ndarray | None = None  # Shape (n_strains,) - accumulated stress / damage

    def __post_init__(self):
        assert self.genomes.ndim == 2
        assert self.genomes.shape[1] == NUM_GENES
        assert len(self.populations) == len(self.genomes)
        if self.lineage_ages is None:
            self.lineage_ages = np.zeros(len(self.populations), dtype=np.float64)
        if self.damage_loads is None:
            self.damage_loads = np.zeros(len(self.populations), dtype=np.float64)
        self.lineage_ages = np.asarray(self.lineage_ages, dtype=np.float64)
        self.damage_loads = np.asarray(self.damage_loads, dtype=np.float64)
        assert len(self.lineage_ages) == len(self.genomes)
        assert len(self.damage_loads) == len(self.genomes)

    @property
    def n_strains(self) -> int:
        return len(self.populations)

    @property
    def total_population(self) -> float:
        return float(np.sum(self.populations))

    def clone(self) -> StrainPopulation:
        return StrainPopulation(
            genomes=self.genomes.copy(),
            populations=self.populations.copy(),
            lineage_ages=self.lineage_ages.copy(),
            damage_loads=self.damage_loads.copy(),
        )

    @classmethod
    def create_initial(
        cls,
        resistant_fraction: float = 0.0,
        initial_population: float = 1e6,
        n_susceptible_strains: int = 3,
        n_resistant_strains: int = 2,
        rng: np.random.Generator = None,
    ) -> StrainPopulation:
        """Create initial population with optional resistance."""
        if rng is None:
            rng = np.random.default_rng()

        strains = []
        pops = []

        # Susceptible strains
        sus_pop = initial_population * (1 - resistant_fraction)
        for i in range(n_susceptible_strains):
            genome = create_wild_type_genome()
            # Add small variation
            genome += rng.normal(0, 0.02, NUM_GENES).astype(np.float32)
            genome = np.clip(genome, 0.0, 1.0)
            strains.append(genome)
            pops.append(sus_pop / n_susceptible_strains)

        # Resistant strains (if any)
        if resistant_fraction > 0 and n_resistant_strains > 0:
            res_pop = initial_population * resistant_fraction
            for i in range(n_resistant_strains):
                resistance_level = 0.3 + rng.random() * 0.4
                genome = create_resistant_genome(resistance_level)
                genome += rng.normal(0, 0.02, NUM_GENES).astype(np.float32)
                genome = np.clip(genome, 0.0, 1.0)
                strains.append(genome)
                pops.append(res_pop / n_resistant_strains)

        genomes = np.array(strains, dtype=np.float32)
        populations = np.array(pops, dtype=np.float64)
        lineage_ages = np.zeros(len(populations), dtype=np.float64)
        damage_loads = np.zeros(len(populations), dtype=np.float64)

        return cls(
            genomes=genomes,
            populations=populations,
            lineage_ages=lineage_ages,
            damage_loads=damage_loads,
        )


def mutate_population(
    population: StrainPopulation,
    config: SimulationConfig,
    abx_stress: float,
    rng: np.random.Generator,
) -> StrainPopulation:
    """
    Apply mutations to population, potentially creating new strains.

    Args:
        population: Current population
        config: Simulation config
        abx_stress: Antibiotic stress level (0-1), increases mutation rate
        rng: Random generator

    Returns:
        Updated population (may have new strains)
    """
    genomes = population.genomes.copy()
    populations = population.populations.copy()
    lineage_ages = population.lineage_ages.copy()
    damage_loads = population.damage_loads.copy()

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

            # Transfer small fraction to mutant
            transfer_fraction = 0.01 * n_mutations
            transfer_pop = populations[i] * transfer_fraction
            populations[i] -= transfer_pop

            new_strains.append(mutant)
            new_pops.append(transfer_pop)
            new_ages.append(max(0.0, lineage_ages[i] * 0.5))
            new_damage.append(max(0.0, damage_loads[i] * 0.7))

    # Add new strains
    if new_strains:
        genomes = np.vstack([genomes, np.array(new_strains, dtype=np.float32)])
        populations = np.concatenate([populations, np.array(new_pops)])
        lineage_ages = np.concatenate([lineage_ages, np.array(new_ages, dtype=np.float64)])
        damage_loads = np.concatenate([damage_loads, np.array(new_damage, dtype=np.float64)])

    return StrainPopulation(
        genomes=genomes,
        populations=populations,
        lineage_ages=lineage_ages,
        damage_loads=damage_loads,
    )


def horizontal_gene_transfer(
    population: StrainPopulation, config: SimulationConfig, rng: np.random.Generator
) -> StrainPopulation:
    """
    Simulate horizontal gene transfer between strains.

    Primarily transfers resistance and persistence genes.
    """
    if population.n_strains < 2:
        return population

    genomes = population.genomes.copy()
    populations = population.populations.copy()
    lineage_ages = population.lineage_ages.copy()
    damage_loads = population.damage_loads.copy()

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

    if new_strains:
        genomes = np.vstack([genomes, np.array(new_strains, dtype=np.float32)])
        populations = np.concatenate([populations, np.array(new_pops)])
        lineage_ages = np.concatenate([lineage_ages, np.array(new_ages, dtype=np.float64)])
        damage_loads = np.concatenate([damage_loads, np.array(new_damage, dtype=np.float64)])

    return StrainPopulation(
        genomes=genomes,
        populations=populations,
        lineage_ages=lineage_ages,
        damage_loads=damage_loads,
    )


def apply_demographic_stochasticity(
    populations: np.ndarray,
    config: SimulationConfig,
    rng: np.random.Generator,
) -> np.ndarray:
    """Realize finite-population noise so small strains can die out by chance."""
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
    immune_status: str,
    rng: np.random.Generator,
) -> StrainPopulation:
    """
    Apply selection based on fitness.

    Grows/shrinks strain populations based on relative fitness.
    """
    if population.n_strains == 0:
        return population

    genomes = population.genomes
    populations = population.populations.copy()
    lineage_ages = population.lineage_ages.copy()
    damage_loads = population.damage_loads.copy()

    # Compute fitness for all strains
    fitness = compute_fitness(
        genomes,
        abx_class=abx_class,
        dose_level=dose_level,
        adherence=adherence,
        immune_strength=immune_strength,
        immune_status=immune_status,
    )

    # Growth based on fitness
    # Relative fitness affects growth rate
    mean_fitness = np.mean(fitness) if len(fitness) > 0 else 0.5
    relative_fitness = fitness / (mean_fitness + 1e-6)

    # Apply selection strength
    selection_factor = np.power(relative_fitness, config.selection_strength)

    # Population dynamics
    growth = config.growth_rate_per_step * selection_factor * fitness
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

    environmental_stress = 1.0 - fitness
    active_dormancy = dormancy_capacity * (0.2 + 0.8 * environmental_stress)

    growth *= np.clip(1.0 - config.dormancy_growth_penalty * active_dormancy, 0.2, 1.0)

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

    death = config.death_rate_per_step * (1.0 / (fitness + 0.1)) + turnover_pressure

    # Net growth
    net_growth = growth - death
    populations = np.maximum(populations * (1.0 + net_growth), 0.0)
    populations = apply_demographic_stochasticity(populations, config, rng)

    # Apply carrying capacity (logistic)
    total = np.sum(populations)
    if total > config.carrying_capacity:
        populations *= config.carrying_capacity / total

    return StrainPopulation(
        genomes=genomes,
        populations=populations,
        lineage_ages=lineage_ages,
        damage_loads=damage_loads,
    )


def consolidate_strains(population: StrainPopulation, config: SimulationConfig) -> StrainPopulation:
    """
    Remove very small strains and merge similar ones to limit strain count.
    """
    # Remove strains below threshold
    mask = population.populations > config.strain_prune_threshold

    if not np.any(mask):
        return StrainPopulation(
            genomes=np.empty((0, NUM_GENES), dtype=np.float32),
            populations=np.empty(0, dtype=np.float64),
            lineage_ages=np.empty(0, dtype=np.float64),
            damage_loads=np.empty(0, dtype=np.float64),
        )

    genomes = population.genomes[mask]
    populations = population.populations[mask]
    lineage_ages = population.lineage_ages[mask]
    damage_loads = population.damage_loads[mask]

    # If too many strains, keep the largest ones
    if len(populations) > config.max_strains:
        indices = np.argsort(populations)[-config.max_strains :]
        genomes = genomes[indices]
        populations = populations[indices]
        lineage_ages = lineage_ages[indices]
        damage_loads = damage_loads[indices]

    return StrainPopulation(
        genomes=genomes,
        populations=populations,
        lineage_ages=lineage_ages,
        damage_loads=damage_loads,
    )


def simulate_day(
    population: StrainPopulation,
    abx_class: str,
    dose_level: str,
    adherence: float,
    immune_strength: float,
    immune_status: str,
    config: SimulationConfig = None,
    seed: int = None,
) -> Tuple[StrainPopulation, Dict[str, Any]]:
    """
    Run 12 simulation steps for one day.

    Args:
        population: Starting population
        abx_class: Antibiotic class being administered
        dose_level: "low", "std", or "high"
        adherence: Patient adherence (0-1)
        immune_strength: Patient immune strength
        immune_status: "normal" or "suppressed"
        config: Simulation configuration
        seed: Random seed for reproducibility

    Returns:
        Tuple of (final_population, step_history)
    """
    if config is None:
        config = SimulationConfig()

    rng = np.random.default_rng(seed)

    # Compute ABX stress level
    profile = ABX_PROFILES.get(abx_class, ABX_PROFILES["none"])
    abx_stress = profile.base_kill_rate * adherence if profile.base_kill_rate > 0 else 0.0

    pop = population.clone()
    history = []

    for step in range(config.steps_per_day):
        # 1. Selection (growth/death based on fitness)
        pop = selection_step(
            pop, config, abx_class, dose_level, adherence, immune_strength, immune_status, rng
        )

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
    immune_status: str,
    config: SimulationConfig = None,
) -> float:
    """
    Compute probability of bacterial clearance (C -> S transition).

    Based on population size and immune factors.
    """
    if config is None:
        config = SimulationConfig()

    total_pop = population.total_population

    if total_pop <= 0:
        return 1.0

    if total_pop < config.clearance_threshold:
        return 0.95  # Very likely to clear

    if total_pop < config.min_population:
        base_prob = 0.3
    else:
        # Logistic function: low clearance at high pop
        ratio = total_pop / config.carrying_capacity
        base_prob = 0.02 * (1.0 - ratio)

    # Immune modulation
    immune_mult = immune_strength
    if immune_status == "suppressed":
        immune_mult *= 0.3

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


def get_dominant_strain(population: StrainPopulation) -> Tuple[np.ndarray, str]:
    """Get the genome and genotype of the dominant strain."""
    if population.n_strains == 0:
        genome = create_wild_type_genome()
        return genome, "S"

    idx = np.argmax(population.populations)
    genome = population.genomes[idx]
    genotype = classify_genotype(genome)

    return genome, genotype


def population_to_response(
    population: StrainPopulation,
    immune_strength: float,
    immune_status: str,
    config: SimulationConfig = None,
) -> Dict[str, Any]:
    """
    Convert final population state to micro response format.

    Returns dict compatible with Patient.apply_micro_response()
    """
    if config is None:
        config = SimulationConfig()

    dominant_genome, dominant_genotype = get_dominant_strain(population)
    resistant_fraction = compute_resistant_fraction(population.genomes, population.populations)

    # Compute derived effects from dominant strain
    transmissibility = float(compute_transmissibility(dominant_genome))
    lethality = float(compute_lethality(dominant_genome))
    severity = float(compute_severity(dominant_genome))
    p_clearance = compute_clearance_probability(population, immune_strength, immune_status, config)

    return {
        "updated_state": {
            "resistant_fraction": resistant_fraction,
            "dominant_genotype": dominant_genotype,
        },
        "derived_effects": {
            "relative_transmissibility": transmissibility,
            "lethality_modifier": lethality,
            "severity_modifier": severity,
            "p_clearance": p_clearance,
        },
        "population_stats": {
            "total_population": population.total_population,
            "n_strains": population.n_strains,
        },
    }
