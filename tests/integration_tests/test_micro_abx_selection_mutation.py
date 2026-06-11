from __future__ import annotations

import numpy as np

from mss.simulation.micro import (
    SimulationConfig,
    StrainPopulation,
    create_resistant_genome,
    create_wild_type_genome,
    mutate_population,
    selection_step,
)


def _two_strain_population() -> StrainPopulation:
    return StrainPopulation(
        genomes=np.vstack(
            [
                create_wild_type_genome(),
                create_resistant_genome(0.8),
            ]
        ).astype(np.float32),
        populations=np.array([1_000_000.0, 1_000_000.0], dtype=np.float64),
        strain_namespace="selection_test",
    )


def test_antibiotic_pressure_amplifies_directional_selection():
    base_config = SimulationConfig(
        selection_strength=1.0,
        abx_selection_pressure_multiplier=0.0,
        stochastic_threshold=0.0,
        stochastic_noise_scale=0.0,
    )
    boosted_config = SimulationConfig(
        selection_strength=1.0,
        abx_selection_pressure_multiplier=3.0,
        stochastic_threshold=0.0,
        stochastic_noise_scale=0.0,
    )

    without_boost = selection_step(
        _two_strain_population(),
        base_config,
        abx_class="beta_lactam",
        dose_level="std",
        adherence=1.0,
        immune_strength=1.0,
        rng=np.random.default_rng(1),
    )
    with_boost = selection_step(
        _two_strain_population(),
        boosted_config,
        abx_class="beta_lactam",
        dose_level="std",
        adherence=1.0,
        immune_strength=1.0,
        rng=np.random.default_rng(1),
    )

    ratio_without_boost = without_boost.populations[1] / without_boost.populations[0]
    ratio_with_boost = with_boost.populations[1] / with_boost.populations[0]

    assert ratio_with_boost > ratio_without_boost


def test_explicit_antibiotic_kill_is_resistance_dependent():
    config = SimulationConfig(
        selection_strength=0.0,
        abx_selection_pressure_multiplier=0.0,
        antibiotic_kill_scale=0.08,
        growth_rate_per_step=0.0,
        base_damage_per_step=0.0,
        replication_damage_factor=0.0,
        stress_damage_factor=0.0,
        repair_rate_per_step=0.0,
        age_mortality_scale=0.0,
        damage_mortality_scale=0.0,
        stochastic_threshold=0.0,
        stochastic_noise_scale=0.0,
    )

    result = selection_step(
        _two_strain_population(),
        config,
        abx_class="beta_lactam",
        dose_level="std",
        adherence=1.0,
        immune_strength=1.0,
        rng=np.random.default_rng(2),
    )

    susceptible_loss = 1_000_000.0 - result.populations[0]
    resistant_loss = 1_000_000.0 - result.populations[1]

    assert susceptible_loss > resistant_loss


def test_mutant_transfer_fraction_controls_founder_population():
    population = StrainPopulation(
        genomes=np.array([create_wild_type_genome()], dtype=np.float32),
        populations=np.array([1_000_000.0], dtype=np.float64),
        strain_namespace="mutation_test",
    )
    low_transfer = SimulationConfig(
        base_mutation_rate=1.0,
        mutation_std=0.05,
        mutant_transfer_fraction=0.01,
    )
    high_transfer = SimulationConfig(
        base_mutation_rate=1.0,
        mutation_std=0.05,
        mutant_transfer_fraction=0.03,
    )

    low_result = mutate_population(
        population.clone(),
        low_transfer,
        abx_stress=0.0,
        rng=np.random.default_rng(3),
    )
    high_result = mutate_population(
        population.clone(),
        high_transfer,
        abx_stress=0.0,
        rng=np.random.default_rng(3),
    )

    assert low_result.n_strains > 1
    assert high_result.n_strains == low_result.n_strains
    assert np.sum(high_result.populations[1:]) > np.sum(low_result.populations[1:])
