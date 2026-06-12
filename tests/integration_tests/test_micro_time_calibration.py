from __future__ import annotations

import numpy as np
import pandas as pd

from tests import micro_config

from mss.simulation.micro import (
    MicroCalibrationScenario,
    TimeScaleChange,
    compound_fraction_per_step_for_new_resolution,
    poisson_intensity_per_step_for_new_resolution,
    probability_per_step_for_new_resolution,
    rescale_micro_config_for_step_duration,
    run_micro_time_scale_ensemble,
    summarize_ensemble,
)


def test_probability_rescaling_preserves_daily_event_probability():
    old_p = 0.06
    new_p = probability_per_step_for_new_resolution(
        old_p,
        reference_opportunities_per_day=12,
        target_opportunities_per_day=24,
    )

    old_daily = 1.0 - (1.0 - old_p) ** 12
    new_daily = 1.0 - (1.0 - new_p) ** 24

    assert np.isclose(new_daily, old_daily)
    assert new_p < old_p


def test_poisson_intensity_rescaling_preserves_daily_expectation():
    change = TimeScaleChange(reference_steps_per_day=12, target_steps_per_day=24)

    old_rate = 0.012
    new_rate = poisson_intensity_per_step_for_new_resolution(old_rate, change)

    assert np.isclose(old_rate * 12, new_rate * 24)


def test_compound_fraction_rescaling_preserves_daily_factor():
    change = TimeScaleChange(reference_steps_per_day=12, target_steps_per_day=24)

    old_growth = 0.18
    new_growth = compound_fraction_per_step_for_new_resolution(old_growth, change)

    assert np.isclose((1.0 + old_growth) ** 12, (1.0 + new_growth) ** 24)
    assert new_growth < old_growth


def test_rescale_micro_config_to_alternative_resolution_changes_step_sensitive_values():
    config = micro_config(steps_per_day=12)
    rescaled = rescale_micro_config_for_step_duration(config, target_steps_per_day=24)

    assert rescaled.steps_per_day == 24
    assert rescaled.base_mutation_rate == config.base_mutation_rate / 2
    assert rescaled.lifecycle_half_life_steps == config.lifecycle_half_life_steps * 2
    assert rescaled.carrying_capacity == config.carrying_capacity
    assert rescaled.selection_strength == config.selection_strength


def test_micro_time_scale_ensemble_summary_has_flat_metric_columns():
    config = micro_config(
        steps_per_day=2,
        max_strains=10,
        base_mutation_rate=0.0,
        base_hgt_rate=0.0,
        strain_prune_threshold=1.0,
    )
    scenario = MicroCalibrationScenario(n_days=1, n_seeds=1, resistant_fraction=0.2)

    rows = run_micro_time_scale_ensemble(config, scenario, label="test")
    summary = summarize_ensemble(rows)

    assert isinstance(rows, pd.DataFrame)
    assert rows.loc[0, "active_window_hours"] == 12.0
    assert rows.loc[0, "step_duration_hours"] == 6.0
    assert list(summary["label"]) == ["test"]
    assert "total_population_mean" in summary.columns
    assert "p_clearance_std" in summary.columns
