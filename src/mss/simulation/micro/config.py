"""Strict configuration models and YAML helpers for the micro simulation."""

from __future__ import annotations

from dataclasses import dataclass, fields
from typing import Any, Mapping


@dataclass
class SimulationConfig:
    """Complete biological configuration for the within-host micro simulation.

    The class intentionally has no field defaults. Production and calibration
    runs must provide every field explicitly through YAML so parameter values are
    visible in the run configuration rather than hidden in Python code.

    Attributes:
        steps_per_day: Number of micro steps simulated per macro day.
        max_strains: Maximum number of active strains retained after pruning.
        carrying_capacity: Logistic carrying capacity for the within-host population.
        min_population: Extinction floor below which carriage is cleared.
        clearance_threshold: Population threshold for near-certain clearance.
        base_mutation_rate: Baseline mutation intensity per gene per step.
        mutation_std: Standard deviation of mutation effects on gene values.
        stress_mutation_boost: Multiplier applied to mutation intensity under ABX stress.
        mutant_transfer_fraction: Parent-population fraction assigned to new mutants.
        base_hgt_rate: Baseline horizontal gene transfer probability per step.
        hgt_gene_transfer_prob: Probability that a transferable gene is blended in an HGT event.
        selection_strength: Exponent applied to relative strain fitness.
        abx_selection_pressure_multiplier: ABX-dependent multiplier for selection strength.
        antibiotic_kill_scale: Per-step explicit antibiotic kill scale.
        growth_rate_per_step: Baseline per-step growth rate before modifiers.
        death_rate_per_step: Baseline per-step death rate.
        death_fitness_floor: Fitness floor used in the baseline death calculation.
        strain_prune_threshold: Absolute population threshold for strain pruning.
        base_damage_per_step: Background damage accumulated per step.
        replication_damage_factor: Damage contribution from replication pressure.
        stress_damage_factor: Damage contribution from environmental stress.
        repair_rate_per_step: Per-step damage repair rate.
        age_mortality_scale: Mortality contribution from lineage age.
        damage_mortality_scale: Mortality contribution from accumulated damage.
        lifecycle_half_life_steps: Scale converting lineage age into turnover pressure.
        max_damage_load: Saturation point for damage-load effects.
        dormancy_growth_penalty: Growth penalty applied to dormant/persister-like cells.
        synergy_repair_dormancy_bonus: Repair bonus from combined repair and dormancy traits.
        synergy_stress_tolerance_bonus: Tolerance bonus from stress response and damage tolerance.
        stochastic_threshold: Population threshold below which Poisson demographic noise is used.
        stochastic_noise_scale: Normal-noise scale for large populations.
        founder_pool_size: Number of reusable founder strains to generate.
        founder_pool_seed: First deterministic seed used for founder generation.
        founder_pool_gene_noise_std: Gene noise applied when generating founder strains.
        gene_presence_threshold: Threshold for binary gene-presence summaries.
    """

    steps_per_day: int
    max_strains: int
    carrying_capacity: float
    min_population: float
    clearance_threshold: float
    base_mutation_rate: float
    mutation_std: float
    stress_mutation_boost: float
    mutant_transfer_fraction: float
    base_hgt_rate: float
    hgt_gene_transfer_prob: float
    selection_strength: float
    abx_selection_pressure_multiplier: float
    antibiotic_kill_scale: float
    growth_rate_per_step: float
    death_rate_per_step: float
    death_fitness_floor: float
    strain_prune_threshold: float
    base_damage_per_step: float
    replication_damage_factor: float
    stress_damage_factor: float
    repair_rate_per_step: float
    age_mortality_scale: float
    damage_mortality_scale: float
    lifecycle_half_life_steps: float
    max_damage_load: float
    dormancy_growth_penalty: float
    synergy_repair_dormancy_bonus: float
    synergy_stress_tolerance_bonus: float
    stochastic_threshold: float
    stochastic_noise_scale: float
    founder_pool_size: int
    founder_pool_seed: int
    founder_pool_gene_noise_std: float
    gene_presence_threshold: float


@dataclass(frozen=True)
class MicroRuntimeConfig:
    """Technical runner settings for the micro simulation.

    Attributes:
        workers: Number of worker processes to use, or None for CPU-count defaulting.
    """

    workers: int | None


def simulation_config_fields() -> set[str]:
    """Return the required biological micro configuration keys.

    Returns:
        Set of field names that must be present in a strict micro YAML block.
    """
    return {field.name for field in fields(SimulationConfig)}


def build_micro_config(raw: Mapping[str, Any], *, source: str = "micro") -> SimulationConfig:
    """Build a complete micro configuration from a YAML mapping.

    The technical ``workers`` key is ignored here because it controls execution,
    not biology. Every ``SimulationConfig`` field must be present and no unknown
    biological keys are accepted.

    Args:
        raw: YAML mapping containing the micro configuration block.
        source: Human-readable source label used in validation errors.

    Returns:
        Complete typed simulation configuration.

    Raises:
        TypeError: If ``raw`` is not a mapping.
        ValueError: If required keys are missing or unknown keys are present.
    """
    if not isinstance(raw, Mapping):
        raise TypeError(f"{source} must be a mapping.")

    values = dict(raw)
    values.pop("workers", None)
    expected = simulation_config_fields()
    actual = set(values)
    missing = sorted(expected - actual)
    unknown = sorted(actual - expected)

    if missing or unknown:
        parts: list[str] = []
        if missing:
            parts.append(f"missing keys: {', '.join(missing)}")
        if unknown:
            parts.append(f"unknown keys: {', '.join(unknown)}")
        raise ValueError(f"Incomplete {source} config ({'; '.join(parts)}).")

    return SimulationConfig(**values)


def parse_micro_runtime_config(
    raw: Mapping[str, Any], *, source: str = "micro"
) -> MicroRuntimeConfig:
    """Extract technical micro-runner settings from a YAML micro block.

    Args:
        raw: YAML mapping containing the micro configuration block.
        source: Human-readable source label used in validation errors.

    Returns:
        Parsed runtime-only configuration.

    Raises:
        TypeError: If ``raw`` is not a mapping.
    """
    if not isinstance(raw, Mapping):
        raise TypeError(f"{source} must be a mapping.")
    workers = raw.get("workers")
    return MicroRuntimeConfig(workers=None if workers is None else int(workers))
