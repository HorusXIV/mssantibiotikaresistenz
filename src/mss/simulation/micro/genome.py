"""
Bacterial genome representation and fitness calculation.

Uses NumPy arrays for vectorized operations across strain populations.
Gene values are normalized floats in [0.0, 1.0].
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Dict

import numpy as np


class GeneIndex(int, Enum):
    """Indices for genome array access."""

    # Metabolism
    GROWTH_BASE = 0
    METABOLIC_OPTIMIZATION = 1
    # Resistance mechanisms
    EFFLUX_PUMPS = 2
    TARGET_MODIFICATION = 3
    PERMEABILITY_REDUCTION = 4
    # Virulence / Survival
    VIRULENCE = 5
    STEALTH = 6
    ADHESION = 7
    # Evolvability
    MUTATION_RATE_MODIFIER = 8
    HGT_COMPETENCE = 9
    # Lifecycle / persistence
    DNA_REPAIR = 10
    DORMANCY_PROPENSITY = 11
    STRESS_RESPONSE = 12
    DAMAGE_TOLERANCE = 13


NUM_GENES = len(GeneIndex)


@dataclass
class ResistanceCosts:
    """Base fitness costs for resistance genes."""

    efflux_pumps: float = 0.15
    target_modification: float = 0.12
    permeability_reduction: float = 0.08

    def as_array(self) -> np.ndarray:
        """Return costs as array aligned with GeneIndex."""
        costs = np.zeros(NUM_GENES, dtype=np.float32)
        costs[GeneIndex.EFFLUX_PUMPS] = self.efflux_pumps
        costs[GeneIndex.TARGET_MODIFICATION] = self.target_modification
        costs[GeneIndex.PERMEABILITY_REDUCTION] = self.permeability_reduction
        return costs


@dataclass
class ABXProfile:
    """Antibiotic class effectiveness against resistance mechanisms."""

    name: str
    # How much each resistance mechanism protects (0-1 multiplier on survival)
    efflux_efficacy: float = 0.5  # How well efflux pumps work against this ABX
    target_mod_efficacy: float = 0.5  # How well target modification works
    permeability_efficacy: float = 0.3  # How well reduced permeability works
    base_kill_rate: float = 0.8  # Base killing rate at standard dose


# Predefined antibiotic profiles (illustrative model values; efficacy/kill in 0-1).
ABX_PROFILES: Dict[str, ABXProfile] = {
    "none": ABXProfile("none", 0.0, 0.0, 0.0, 0.0),
    # Beta-lactams struggle against efflux but are destroyed by modification
    # (or bypass via PBP2a). Permeability changes hit them hard.
    "beta_lactam": ABXProfile("beta_lactam", 0.1, 0.9, 0.6, 0.75),
    # Fluoroquinolones are largely unaffected by porins, but target mod
    # (gyrase mutations) is huge. Efflux is moderate.
    "fluoroquinolone": ABXProfile("fluoroquinolone", 0.5, 0.8, 0.1, 0.80),
    # Aminoglycosides are big targets for modifying enzymes and permeability
    # changes, but not classic efflux pumps.
    "aminoglycoside": ABXProfile("aminoglycoside", 0.2, 0.7, 0.7, 0.70),
    # Macrolides are heavily effluxed and subject to ribosomal modification.
    "macrolide": ABXProfile("macrolide", 0.7, 0.7, 0.2, 0.65),
    # Tetracyclines are the poster child for efflux pumps.
    "tetracycline": ABXProfile("tetracycline", 0.9, 0.2, 0.2, 0.60),
    # Glycopeptides are too large to be pumped out (Efflux = 0),
    # but modification (D-ala-D-lac) makes them completely useless.
    "glycopeptide": ABXProfile("glycopeptide", 0.0, 0.95, 0.4, 0.85),
}

DOSE_MULTIPLIERS = {
    "low": 0.6,
    "std": 1.0,
    "high": 1.4,
}


def create_wild_type_genome() -> np.ndarray:
    """Create a susceptible wild-type genome."""
    genome = np.zeros(NUM_GENES, dtype=np.float32)
    genome[GeneIndex.GROWTH_BASE] = 0.8  # Good baseline growth
    genome[GeneIndex.METABOLIC_OPTIMIZATION] = 0.1
    genome[GeneIndex.VIRULENCE] = 0.3
    genome[GeneIndex.STEALTH] = 0.2
    genome[GeneIndex.ADHESION] = 0.4
    genome[GeneIndex.MUTATION_RATE_MODIFIER] = 0.1
    genome[GeneIndex.HGT_COMPETENCE] = 0.2
    genome[GeneIndex.DNA_REPAIR] = 0.2
    genome[GeneIndex.DORMANCY_PROPENSITY] = 0.1
    genome[GeneIndex.STRESS_RESPONSE] = 0.2
    genome[GeneIndex.DAMAGE_TOLERANCE] = 0.15
    # Resistance genes start near zero
    genome[GeneIndex.EFFLUX_PUMPS] = 0.05
    genome[GeneIndex.TARGET_MODIFICATION] = 0.05
    genome[GeneIndex.PERMEABILITY_REDUCTION] = 0.05
    return genome


def create_resistant_genome(resistance_level: float = 0.5) -> np.ndarray:
    """Create a partially resistant genome."""
    genome = create_wild_type_genome()
    genome[GeneIndex.EFFLUX_PUMPS] = resistance_level * 0.8
    genome[GeneIndex.TARGET_MODIFICATION] = resistance_level * 0.6
    genome[GeneIndex.PERMEABILITY_REDUCTION] = resistance_level * 0.4
    genome[GeneIndex.METABOLIC_OPTIMIZATION] = resistance_level * 0.5  # Compensatory
    genome[GeneIndex.DNA_REPAIR] = 0.15 + resistance_level * 0.35
    genome[GeneIndex.STRESS_RESPONSE] = 0.20 + resistance_level * 0.40
    genome[GeneIndex.DAMAGE_TOLERANCE] = 0.10 + resistance_level * 0.35
    genome[GeneIndex.DORMANCY_PROPENSITY] = 0.05 + resistance_level * 0.20
    genome[GeneIndex.GROWTH_BASE] = 0.7  # Slightly reduced
    return genome


def compute_resistance_costs(genomes: np.ndarray, costs: ResistanceCosts = None) -> np.ndarray:
    """
    Compute net fitness cost from resistance genes.

    Args:
        genomes: Shape (n_strains, NUM_GENES) or (NUM_GENES,)
        costs: ResistanceCosts parameters

    Returns:
        Net costs array of shape (n_strains,) or scalar
    """
    if costs is None:
        costs = ResistanceCosts()

    single = genomes.ndim == 1
    if single:
        genomes = genomes.reshape(1, -1)

    cost_array = costs.as_array()

    # Raw cost = sum(gene_value * gene_cost)
    raw_costs = np.sum(genomes * cost_array, axis=1)

    # Metabolic optimization reduces costs (up to 80% reduction at max)
    optimization = genomes[:, GeneIndex.METABOLIC_OPTIMIZATION]
    cost_reduction = optimization * 0.8  # Max 80% cost reduction

    net_costs = raw_costs * (1.0 - cost_reduction)

    return net_costs[0] if single else net_costs


def compute_abx_survival(
    genomes: np.ndarray, abx_class: str, dose_level: str, adherence: float
) -> np.ndarray:
    """
    Compute survival factor under antibiotic pressure.

    Args:
        genomes: Shape (n_strains, NUM_GENES) or (NUM_GENES,)
        abx_class: Antibiotic class name
        dose_level: "low", "std", or "high"
        adherence: Patient adherence 0-1

    Returns:
        Survival factor array (0-1), shape (n_strains,) or scalar
    """
    profile = ABX_PROFILES.get(abx_class, ABX_PROFILES["none"])

    if profile.base_kill_rate == 0.0:
        # No antibiotic pressure
        single = genomes.ndim == 1
        n = 1 if single else genomes.shape[0]
        return 1.0 if single else np.ones(n, dtype=np.float32)

    single = genomes.ndim == 1
    if single:
        genomes = genomes.reshape(1, -1)

    dose_mult = DOSE_MULTIPLIERS.get(dose_level, 1.0)
    effective_kill = profile.base_kill_rate * dose_mult * adherence

    # Resistance protection
    efflux = genomes[:, GeneIndex.EFFLUX_PUMPS]
    target_mod = genomes[:, GeneIndex.TARGET_MODIFICATION]
    permeability = genomes[:, GeneIndex.PERMEABILITY_REDUCTION]

    # Combined resistance effect (multiplicative protection)
    weighted_resistance = (
        efflux * profile.efflux_efficacy
        + target_mod * profile.target_mod_efficacy
        + permeability * profile.permeability_efficacy
    )
    protection = weighted_resistance / 1
    # Normalize and cap at 0.99
    protection = np.clip(protection, 0.0, 0.99)

    # Survival = 1 - kill_rate * (1 - protection) ** 2
    susceptibility = (1.0 - protection) ** 2
    abx_survival = 1.0 - effective_kill * susceptibility
    survival = np.clip(abx_survival, 0.01, 1.0)

    return survival[0] if single else survival


def compute_immune_survival(genomes: np.ndarray, immune_strength: float) -> np.ndarray:
    """
    Compute survival factor against immune system.

    Args:
        genomes: Shape (n_strains, NUM_GENES) or (NUM_GENES,)
        immune_strength: Patient immune strength multiplier (lower = weaker
            immune competence; immunocompromised patients carry a low value)

    Returns:
        Survival factor array (0-1)
    """
    single = genomes.ndim == 1
    if single:
        genomes = genomes.reshape(1, -1)

    # Base immune clearance rate
    base_clearance = 0.15 * immune_strength

    # Stealth gene reduces immune detection
    stealth = genomes[:, GeneIndex.STEALTH]
    evasion = stealth * 0.7  # Up to 70% immune evasion

    survival = 1.0 - base_clearance * (1.0 - evasion)
    survival = np.clip(survival, 0.05, 1.0)

    return survival[0] if single else survival


def compute_fitness(
    genomes: np.ndarray,
    abx_class: str = "none",
    dose_level: str = "std",
    adherence: float = 1.0,
    immune_strength: float = 1.0,
    costs: ResistanceCosts = None,
) -> np.ndarray:
    """
    Compute overall fitness for bacterial genomes.

    Fitness = (growth_base - net_costs) * ABX_survival * Immune_survival

    Args:
        genomes: Shape (n_strains, NUM_GENES) or (NUM_GENES,)
        abx_class: Antibiotic class
        dose_level: Dose level
        adherence: Patient adherence
        immune_strength: Patient immune strength
        costs: Resistance cost parameters

    Returns:
        Fitness array, shape (n_strains,) or scalar
    """
    single = genomes.ndim == 1
    if single:
        genomes = genomes.reshape(1, -1)

    # Base growth
    growth = genomes[:, GeneIndex.GROWTH_BASE]

    # Subtract resistance costs
    net_costs = compute_resistance_costs(genomes, costs)
    base_fitness = np.clip(growth - net_costs, 0.01, 1.0)

    # Multiply by survival factors
    abx_survival = compute_abx_survival(genomes, abx_class, dose_level, adherence)
    immune_survival = compute_immune_survival(genomes, immune_strength)

    fitness = base_fitness * abx_survival * immune_survival
    fitness = np.clip(fitness, 0.001, 1.0)

    return fitness[0] if single else fitness


def compute_transmissibility(genomes: np.ndarray) -> np.ndarray:
    """
    Compute relative transmissibility based on adhesion and virulence.

    Args:
        genomes: Shape (n_strains, NUM_GENES) or (NUM_GENES,)

    Returns:
        Transmissibility multiplier (centered around 1.0)
    """
    single = genomes.ndim == 1
    if single:
        genomes = genomes.reshape(1, -1)

    adhesion = genomes[:, GeneIndex.ADHESION]
    virulence = genomes[:, GeneIndex.VIRULENCE]

    # Adhesion is primary driver, virulence secondary
    transmissibility = 0.5 + adhesion * 0.8 + virulence * 0.2

    return transmissibility[0] if single else transmissibility


def compute_severity(genomes: np.ndarray) -> np.ndarray:
    """
    Compute severity modifier based on virulence and adhesion.

    Args:
        genomes: Shape (n_strains, NUM_GENES) or (NUM_GENES,)

    Returns:
        Severity multiplier (centered around 1.0)
    """
    single = genomes.ndim == 1
    if single:
        genomes = genomes.reshape(1, -1)

    virulence = genomes[:, GeneIndex.VIRULENCE]
    adhesion = genomes[:, GeneIndex.ADHESION]

    severity = 0.6 + virulence * 0.8 + adhesion * 0.2

    return severity[0] if single else severity


def classify_genotype(genome: np.ndarray) -> str:
    """
    Classify a genome into a genotype category.

    Returns: "S" (susceptible), "R1" (low resistance), "R2" (medium), "R3" (high)
    """
    resistance_score = (
        genome[GeneIndex.EFFLUX_PUMPS]
        + genome[GeneIndex.TARGET_MODIFICATION]
        + genome[GeneIndex.PERMEABILITY_REDUCTION]
    ) / 3.0

    # Model thresholds on the mean resistance score (not clinically calibrated).
    if resistance_score < 0.2:
        return "S"
    elif resistance_score < 0.4:
        return "R1"
    elif resistance_score < 0.7:
        return "R2"
    else:
        return "R3"


def compute_resistant_fraction(genomes: np.ndarray, populations: np.ndarray) -> float:
    """
    Compute fraction of population that is resistant.

    Args:
        genomes: Shape (n_strains, NUM_GENES)
        populations: Shape (n_strains,) - population sizes

    Returns:
        Fraction of total population with resistance score >= 0.3
    """
    resistance_scores = (
        genomes[:, GeneIndex.EFFLUX_PUMPS]
        + genomes[:, GeneIndex.TARGET_MODIFICATION]
        + genomes[:, GeneIndex.PERMEABILITY_REDUCTION]
    ) / 3.0

    total_pop = np.sum(populations)
    if total_pop == 0:
        return 0.0

    resistant_mask = resistance_scores >= 0.3
    resistant_pop = np.sum(populations[resistant_mask])

    return float(resistant_pop / total_pop)
