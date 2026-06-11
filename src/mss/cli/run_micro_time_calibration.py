"""Generate and validate micro time-scale calibration candidates."""

from __future__ import annotations

import argparse
import copy
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from mss.cli.run_coupled_simulation import DEFAULT_CONFIG_PATH, PROJECT_ROOT, load_coupled_settings
from mss.simulation.micro import (
    MicroCalibrationScenario,
    describe_time_scaling,
    rescale_micro_config_for_step_duration,
    run_micro_time_scale_ensemble,
    summarize_ensemble,
)


def _load_raw_config(path: Path) -> dict[str, Any]:
    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict):
        raise TypeError("Top-level config must be a mapping.")
    return raw


def _write_candidate_config(
    raw: dict[str, Any],
    output_path: Path,
    micro_values: dict[str, Any],
) -> None:
    candidate = copy.deepcopy(raw)
    candidate.setdefault("micro", {})
    workers = candidate["micro"].get("workers", None)
    candidate["micro"] = dict(micro_values)
    candidate["micro"]["workers"] = workers
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(yaml.safe_dump(candidate, sort_keys=False), encoding="utf-8")


def _scenario_from_args(args: argparse.Namespace) -> MicroCalibrationScenario:
    return MicroCalibrationScenario(
        n_days=args.n_days,
        n_seeds=args.n_seeds,
        resistant_fraction=args.resistant_fraction,
        dominant_genotype=args.dominant_genotype,
        abx_on=args.abx_on,
        abx_class=args.abx_class,
        dose_level=args.dose_level,
        adherence=args.adherence,
        immune_strength=args.immune_strength,
        initial_population=args.initial_population,
        active_window_hours=args.active_window_hours,
    )


def _print_final_comparison(summary: pd.DataFrame, final_day: int) -> None:
    final = summary[summary["day"] == final_day]
    columns = [
        "label",
        "steps_per_day",
        "total_population_mean",
        "resistant_fraction_mean",
        "p_clearance_mean",
        "n_strains_mean",
    ]
    available = [col for col in columns if col in final.columns]
    print("\nFinal-day ensemble means:")
    print(final[available].to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Validate the 12-step overnight micro window, or rescale micro "
            "parameters if a different active-window resolution is requested."
        )
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG_PATH)
    parser.add_argument("--target-steps-per-day", type=int, default=12)
    parser.add_argument("--reference-steps-per-day", type=int, default=None)
    parser.add_argument("--n-days", type=int, default=30)
    parser.add_argument("--n-seeds", type=int, default=5)
    parser.add_argument("--resistant-fraction", type=float, default=0.9)
    parser.add_argument("--dominant-genotype", type=str, default="R2")
    parser.add_argument("--abx-on", action="store_true")
    parser.add_argument("--abx-class", type=str, default="beta_lactam")
    parser.add_argument("--dose-level", type=str, default="std")
    parser.add_argument("--adherence", type=float, default=0.7)
    parser.add_argument("--immune-strength", type=float, default=0.75)
    parser.add_argument("--initial-population", type=float, default=1e6)
    parser.add_argument("--active-window-hours", type=float, default=12.0)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    settings = load_coupled_settings(args.config)
    raw_config = _load_raw_config(args.config)
    scenario = _scenario_from_args(args)
    reference_config = settings.micro
    target_config = rescale_micro_config_for_step_duration(
        reference_config,
        target_steps_per_day=args.target_steps_per_day,
        reference_steps_per_day=args.reference_steps_per_day,
    )

    output_dir = args.output_dir or (
        PROJECT_ROOT
        / "outputs"
        / (
            datetime.now().strftime("%Y%m%d_%H%M%S")
            + f"_MicroTime_{args.target_steps_per_day}_steps_{args.active_window_hours:g}h"
        )
    )
    data_dir = output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    scaling = describe_time_scaling(reference_config, target_config)
    scaling_path = data_dir / "micro_time_scaling.parquet"
    scaling.to_parquet(scaling_path, index=False)

    candidate_path = data_dir / f"candidate_micro_{args.target_steps_per_day}_steps.yml"
    _write_candidate_config(raw_config, candidate_path, asdict(target_config))

    reference_df = run_micro_time_scale_ensemble(
        reference_config,
        scenario,
        seed_offset=0,
        label="reference",
    )
    target_df = run_micro_time_scale_ensemble(
        target_config,
        scenario,
        seed_offset=0,
        label="target_rescaled",
    )
    ensemble = pd.concat([reference_df, target_df], ignore_index=True)
    summary = summarize_ensemble(ensemble)

    ensemble_path = data_dir / "micro_time_ensemble.parquet"
    summary_path = data_dir / "micro_time_summary.parquet"
    ensemble.to_parquet(ensemble_path, index=False)
    summary.to_parquet(summary_path, index=False)

    print("=" * 72)
    print("Micro time-scale calibration")
    print(f"  Config: {args.config}")
    print(f"  Reference steps/day: {reference_config.steps_per_day}")
    print(f"  Target steps/day: {target_config.steps_per_day}")
    print(f"  Active micro window: {scenario.active_window_hours:g} hours per macro day")
    print(
        f"  Target step duration: {scenario.active_window_hours / target_config.steps_per_day:g} hours"
    )
    print(f"  Scenario days: {scenario.n_days}")
    print(f"  Seeds: {scenario.n_seeds}")
    print(f"  Antibiotics: {'on' if scenario.abx_on else 'off'} ({scenario.abx_class})")
    _print_final_comparison(summary, scenario.n_days)
    print("\nSaved:")
    print(f"  Candidate config -> {candidate_path}")
    print(f"  Scaling table    -> {scaling_path}")
    print(f"  Ensemble rows    -> {ensemble_path}")
    print(f"  Summary          -> {summary_path}")
    print("=" * 72)


if __name__ == "__main__":
    main()
