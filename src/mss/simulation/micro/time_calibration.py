"""Time-scale calibration helpers for the micro simulation.

The micro engine stores several parameters per simulation step. When the
meaning of one step changes, those parameters cannot be copied unchanged:
probabilities and rates must be transformed so the daily process remains
comparable.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from typing import Any

import numpy as np
import pandas as pd

try:
    from .engine import SimulationConfig, StrainPopulation, population_to_response, simulate_day
    from .genome import compute_resistant_fraction
except ImportError:
    from mss.simulation.micro.engine import (
        SimulationConfig,
        StrainPopulation,
        population_to_response,
        simulate_day,
    )
    from mss.simulation.micro.genome import compute_resistant_fraction


@dataclass(frozen=True)
class TimeScaleChange:
    """Mapping between a reference micro resolution and a target resolution."""

    reference_steps_per_day: int
    target_steps_per_day: int

    def __post_init__(self) -> None:
        if self.reference_steps_per_day <= 0:
            raise ValueError("reference_steps_per_day must be positive.")
        if self.target_steps_per_day <= 0:
            raise ValueError("target_steps_per_day must be positive.")

    @property
    def step_scale(self) -> float:
        """Reference-step duration divided by target-step duration."""
        return self.reference_steps_per_day / self.target_steps_per_day


@dataclass(frozen=True)
class MicroCalibrationScenario:
    """Controlled within-host scenario used for time-scale validation."""

    n_days: int = 30
    n_seeds: int = 5
    resistant_fraction: float = 0.9
    dominant_genotype: str = "R2"
    abx_on: bool = False
    abx_class: str = "beta_lactam"
    dose_level: str = "std"
    adherence: float = 0.7
    immune_strength: float = 0.75
    initial_population: float = 1e6
    active_window_hours: float = 12.0

    def __post_init__(self) -> None:
        if self.n_days <= 0:
            raise ValueError("n_days must be positive.")
        if self.n_seeds <= 0:
            raise ValueError("n_seeds must be positive.")
        if not 0.0 <= self.resistant_fraction <= 1.0:
            raise ValueError("resistant_fraction must be in [0, 1].")
        if self.initial_population <= 0.0:
            raise ValueError("initial_population must be positive.")
        if self.active_window_hours <= 0.0:
            raise ValueError("active_window_hours must be positive.")


def probability_per_step_for_new_resolution(
    probability_per_step: float,
    *,
    reference_opportunities_per_day: float,
    target_opportunities_per_day: float,
) -> float:
    """Rescale an event probability while preserving the daily no-event probability."""
    if not 0.0 <= probability_per_step <= 1.0:
        raise ValueError("probability_per_step must be in [0, 1].")
    if reference_opportunities_per_day <= 0.0 or target_opportunities_per_day <= 0.0:
        raise ValueError("opportunities per day must be positive.")

    daily_survival = (1.0 - probability_per_step) ** reference_opportunities_per_day
    return float(1.0 - daily_survival ** (1.0 / target_opportunities_per_day))


def poisson_intensity_per_step_for_new_resolution(
    intensity_per_step: float,
    change: TimeScaleChange,
) -> float:
    """Rescale a Poisson intensity so the expected daily count is unchanged."""
    if intensity_per_step < 0.0:
        raise ValueError("intensity_per_step must be non-negative.")
    return float(intensity_per_step * change.step_scale)


def compound_fraction_per_step_for_new_resolution(
    fraction_per_step: float,
    change: TimeScaleChange,
) -> float:
    """Rescale a fractional process while preserving its daily compound factor."""
    if fraction_per_step <= -1.0:
        raise ValueError("fraction_per_step must be greater than -1.")
    return float(np.expm1(np.log1p(fraction_per_step) * change.step_scale))


def additive_per_step_for_new_resolution(value_per_step: float, change: TimeScaleChange) -> float:
    """Rescale a quantity that is accumulated linearly per step."""
    return float(value_per_step * change.step_scale)


def rescale_micro_config_for_step_duration(
    config: SimulationConfig,
    *,
    target_steps_per_day: int,
    reference_steps_per_day: int | None = None,
) -> SimulationConfig:
    """Return a copy of ``config`` calibrated to a new micro step duration.

    The function treats ``reference_steps_per_day`` as the resolution at which
    the current values were calibrated. If omitted, ``config.steps_per_day`` is
    used. Population thresholds, trait magnitudes, selection shape, and capacity
    limits are left unchanged because they are not per-step rates.
    """
    reference = int(reference_steps_per_day or config.steps_per_day)
    change = TimeScaleChange(
        reference_steps_per_day=reference,
        target_steps_per_day=int(target_steps_per_day),
    )
    hgt_reference_opportunities = max(1.0, reference / 3.0)
    hgt_target_opportunities = max(1.0, target_steps_per_day / 3.0)

    return replace(
        config,
        steps_per_day=change.target_steps_per_day,
        base_mutation_rate=poisson_intensity_per_step_for_new_resolution(
            config.base_mutation_rate,
            change,
        ),
        base_hgt_rate=probability_per_step_for_new_resolution(
            config.base_hgt_rate,
            reference_opportunities_per_day=hgt_reference_opportunities,
            target_opportunities_per_day=hgt_target_opportunities,
        ),
        growth_rate_per_step=compound_fraction_per_step_for_new_resolution(
            config.growth_rate_per_step,
            change,
        ),
        death_rate_per_step=probability_per_step_for_new_resolution(
            config.death_rate_per_step,
            reference_opportunities_per_day=reference,
            target_opportunities_per_day=target_steps_per_day,
        ),
        base_damage_per_step=additive_per_step_for_new_resolution(
            config.base_damage_per_step,
            change,
        ),
        replication_damage_factor=additive_per_step_for_new_resolution(
            config.replication_damage_factor,
            change,
        ),
        stress_damage_factor=additive_per_step_for_new_resolution(
            config.stress_damage_factor,
            change,
        ),
        repair_rate_per_step=additive_per_step_for_new_resolution(
            config.repair_rate_per_step,
            change,
        ),
        age_mortality_scale=additive_per_step_for_new_resolution(
            config.age_mortality_scale,
            change,
        ),
        lifecycle_half_life_steps=float(config.lifecycle_half_life_steps / change.step_scale),
    )


def describe_time_scaling(before: SimulationConfig, after: SimulationConfig) -> pd.DataFrame:
    """Return a table of changed micro parameters and multiplicative factors."""
    rows: list[dict[str, Any]] = []
    before_raw = asdict(before)
    after_raw = asdict(after)
    for name, old_value in before_raw.items():
        new_value = after_raw[name]
        if old_value == new_value:
            continue
        factor = (
            float(new_value) / float(old_value)
            if isinstance(old_value, (int, float)) and float(old_value) != 0.0
            else np.nan
        )
        rows.append(
            {
                "parameter": name,
                "reference_value": old_value,
                "target_value": new_value,
                "factor": factor,
            }
        )
    return pd.DataFrame(
        rows,
        columns=["parameter", "reference_value", "target_value", "factor"],
    )


def run_micro_time_scale_ensemble(
    config: SimulationConfig,
    scenario: MicroCalibrationScenario,
    *,
    seed_offset: int = 0,
    label: str = "candidate",
) -> pd.DataFrame:
    """Run controlled micro episodes and return one row per seed and day."""
    rows: list[dict[str, Any]] = []
    abx_class = scenario.abx_class if scenario.abx_on else "none"

    for seed_index in range(scenario.n_seeds):
        seed = seed_offset + seed_index
        rng = np.random.default_rng(seed)
        population = StrainPopulation.create_initial(
            resistant_fraction=scenario.resistant_fraction,
            dominant_genotype=scenario.dominant_genotype,
            initial_population=scenario.initial_population,
            rng=rng,
            strain_namespace=f"{label}_{seed_index}",
        )

        for day in range(1, scenario.n_days + 1):
            population, _history = simulate_day(
                population=population,
                abx_class=abx_class,
                dose_level=scenario.dose_level,
                adherence=scenario.adherence,
                immune_strength=scenario.immune_strength,
                config=config,
                seed=seed * 100_000 + day,
            )
            response = population_to_response(
                population,
                immune_strength=scenario.immune_strength,
                config=config,
            )
            rows.append(
                {
                    "label": label,
                    "seed": seed,
                    "day": day,
                    "steps_per_day": config.steps_per_day,
                    "active_window_hours": scenario.active_window_hours,
                    "step_duration_hours": scenario.active_window_hours / config.steps_per_day,
                    "total_population": population.total_population,
                    "n_strains": population.n_strains,
                    "resistant_fraction": compute_resistant_fraction(
                        population.genomes,
                        population.populations,
                    ),
                    "p_clearance": response["derived_effects"]["p_clearance"],
                    "relative_transmissibility": response["derived_effects"][
                        "relative_transmissibility"
                    ],
                    "severity_modifier": response["derived_effects"]["severity_modifier"],
                    "dominant_genotype": response["updated_state"]["dominant_genotype"],
                }
            )

    return pd.DataFrame(rows)


def summarize_ensemble(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate daily ensemble output into mean/std calibration diagnostics."""
    metrics = [
        "total_population",
        "n_strains",
        "resistant_fraction",
        "p_clearance",
        "relative_transmissibility",
        "severity_modifier",
    ]
    summary = (
        df.groupby(["label", "steps_per_day", "day"], as_index=False)[metrics]
        .agg(["mean", "std"])
        .reset_index()
    )
    summary.columns = [
        "_".join(str(part) for part in col if part) if isinstance(col, tuple) else str(col)
        for col in summary.columns
    ]
    return summary
