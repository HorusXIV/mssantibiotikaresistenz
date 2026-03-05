"""
Micro-simulation module for within-host bacterial evolution.

This module simulates bacterial population dynamics using an evolutionary
algorithm. Each day is divided into 12 discrete time steps where selection,
mutation, and horizontal gene transfer occur.

Key Components:
- BacterialGenome: Gene representation with 10 normalized attributes
- StrainPopulation: Population model tracking multiple strains
- MicroSimulator: Interface for batch processing patient episodes

Usage:
    from micro_simulation import MicroSimulator, SimulationConfig

    simulator = MicroSimulator()
    response = simulator.process_request(patient.make_micro_request(...))
    patient.apply_micro_response(response)
"""

from .genome import (
    GeneIndex,
    NUM_GENES,
    ResistanceCosts,
    ABXProfile,
    ABX_PROFILES,
    DOSE_MULTIPLIERS,
    create_wild_type_genome,
    create_resistant_genome,
    compute_fitness,
    compute_resistance_costs,
    compute_abx_survival,
    compute_immune_survival,
    compute_transmissibility,
    compute_lethality,
    compute_severity,
    classify_genotype,
    compute_resistant_fraction,
)

from .simulation import (
    SimulationConfig,
    StrainPopulation,
    simulate_day,
    mutate_population,
    horizontal_gene_transfer,
    selection_step,
    compute_clearance_probability,
    get_dominant_strain,
    population_to_response,
)

from .simulator import (
    EpisodeState,
    MicroSimulator,
    run_micro_simulation,
)

__all__ = [
    # Genome
    "GeneIndex",
    "NUM_GENES",
    "ResistanceCosts",
    "ABXProfile",
    "ABX_PROFILES",
    "DOSE_MULTIPLIERS",
    "create_wild_type_genome",
    "create_resistant_genome",
    "compute_fitness",
    "compute_resistance_costs",
    "compute_abx_survival",
    "compute_immune_survival",
    "compute_transmissibility",
    "compute_lethality",
    "compute_severity",
    "classify_genotype",
    "compute_resistant_fraction",
    # Simulation
    "SimulationConfig",
    "StrainPopulation",
    "simulate_day",
    "mutate_population",
    "horizontal_gene_transfer",
    "selection_step",
    "compute_clearance_probability",
    "get_dominant_strain",
    "population_to_response",
    # Simulator interface
    "EpisodeState",
    "MicroSimulator",
    "run_micro_simulation",
]
