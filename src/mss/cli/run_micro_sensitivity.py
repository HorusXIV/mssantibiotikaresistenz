"""Micro-parameter sensitivity analysis (Calibration 7).

For each calibratable YAML parameter in the micro block, sweep its value while
keeping all others fixed at the baseline (from the config), and record how key
within-host metrics respond over time.

Metrics tracked daily (mean ± std over seeds):
  resistant_fraction  -- fraction of within-host population with resistance score >= 0.3
  total_population    -- absolute bacterial count
  p_clearance         -- daily probability of infection clearing
  n_strains           -- number of active strains

Two fixed scenarios per parameter:
  no_abx      -- natural persistence without antibiotic pressure
  beta_lactam -- standard beta-lactam treatment (dose=std, adherence=0.7)

Parameters swept (all non-technical, calibratable micro YAML keys):
  Selection & dynamics:
    selection_strength, death_fitness_floor, growth_rate_per_step,
    death_rate_per_step, strain_prune_threshold
  Mutation & HGT:
    base_mutation_rate, mutation_std, stress_mutation_boost,
    base_hgt_rate, hgt_gene_transfer_prob
  Damage & lifecycle:
    base_damage_per_step, replication_damage_factor, stress_damage_factor,
    repair_rate_per_step, age_mortality_scale, damage_mortality_scale,
    lifecycle_half_life_steps, max_damage_load, dormancy_growth_penalty,
    synergy_repair_dormancy_bonus, synergy_stress_tolerance_bonus
  Population thresholds:
    carrying_capacity, min_population, clearance_threshold

Outputs (outputs/<timestamp>_MicroSensitivity/):
  data/sensitivity_trajectories.parquet  -- daily metrics for all runs
  data/sensitivity_summary.parquet       -- final-day stats per (param, value, scenario)
  figures/sensitivity_<param>.png        -- trajectory plots per parameter
  figures/sensitivity_overview.png       -- sensitivity magnitude overview

Usage:
    uv run python -m mss.cli.run_micro_sensitivity
    uv run python -m mss.cli.run_micro_sensitivity --config config/cal_micro_sensitivity.yml
    uv run python -m mss.cli.run_micro_sensitivity --params selection_strength base_mutation_rate
    uv run python -m mss.cli.run_micro_sensitivity --n-days 90 --n-seeds 5
    uv run python -m mss.cli.run_micro_sensitivity --no-abx-only
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
from datetime import datetime
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.use("Agg")

from mss.cli.run_coupled_simulation import PROJECT_ROOT, load_coupled_settings
from mss.simulation.micro.engine import (
    SimulationConfig,
    StrainPopulation,
    compute_clearance_probability,
    simulate_day,
)
from mss.simulation.micro.genome import compute_resistant_fraction

# ---------------------------------------------------------------------------
# Parameter sweep definitions
# Each entry: (config_field, label, base_value, sweep_values)
# base_value is taken from cal_micro_sensitivity.yml / simulation_realistic_micro.yml
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SweepParam:
    field: str
    label: str
    values: list[float]
    log_scale: bool = False


# fmt: off
_SWEEP_PARAMS: list[SweepParam] = [
    # --- Selection & Population Dynamics ---
    SweepParam("selection_strength",       "Selection strength",          [0.05, 0.15, 0.3, 0.6, 1.2, 2.5, 5.0]),
    SweepParam("death_fitness_floor",      "Death fitness floor",         [0.05, 0.2, 0.5, 1.0, 2.0, 4.0, 8.0]),
    SweepParam("growth_rate_per_step",     "Growth rate / step",          [0.2, 0.5, 0.8, 1.5, 2.5, 4.0, 6.0]),
    SweepParam("death_rate_per_step",      "Death rate / step",           [0.001, 0.003, 0.005, 0.01, 0.025, 0.06, 0.12]),
    SweepParam("strain_prune_threshold",   "Strain prune threshold",      [5.0, 20.0, 50.0, 200.0, 1000.0, 5000.0, 20000.0], log_scale=True),

    # --- Mutation & HGT ---
    SweepParam("base_mutation_rate",       "Base mutation rate",          [0.0005, 0.001, 0.003, 0.006, 0.015, 0.04, 0.10]),
    SweepParam("mutation_std",             "Mutation std",                [0.002, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20]),
    SweepParam("stress_mutation_boost",    "Stress mutation boost",       [2.0, 5.0, 15.0, 30.0, 60.0, 120.0, 250.0]),
    SweepParam("base_hgt_rate",            "HGT base rate",               [0.001, 0.004, 0.008, 0.02, 0.05, 0.12, 0.25]),
    SweepParam("hgt_gene_transfer_prob",   "HGT gene transfer prob",      [0.02, 0.05, 0.12, 0.25, 0.50, 0.70, 0.90]),

    # --- Damage & Lifecycle ---
    SweepParam("base_damage_per_step",     "Base damage / step",          [0.0005, 0.001, 0.002, 0.004, 0.012, 0.030, 0.08]),
    SweepParam("replication_damage_factor","Replication damage factor",   [0.002, 0.007, 0.015, 0.03, 0.06, 0.12, 0.25]),
    SweepParam("stress_damage_factor",     "Stress damage factor",        [0.005, 0.015, 0.03, 0.06, 0.12, 0.25, 0.50]),
    SweepParam("repair_rate_per_step",     "Repair rate / step",          [0.005, 0.02, 0.04, 0.08, 0.15, 0.30, 0.60]),
    SweepParam("age_mortality_scale",      "Age mortality scale",         [0.00005, 0.0002, 0.0005, 0.001, 0.005, 0.02, 0.08]),
    SweepParam("damage_mortality_scale",   "Damage mortality scale",      [0.002, 0.007, 0.012, 0.025, 0.06, 0.15, 0.35]),
    SweepParam("lifecycle_half_life_steps","Lifecycle half-life (steps)", [20.0, 50.0, 100.0, 200.0, 500.0, 2000.0, 8000.0]),
    SweepParam("max_damage_load",          "Max damage load",             [0.5, 1.5, 2.5, 5.0, 10.0, 20.0, 50.0]),
    SweepParam("dormancy_growth_penalty",  "Dormancy growth penalty",     [0.05, 0.15, 0.30, 0.55, 0.75, 0.90, 0.99]),
    SweepParam("synergy_repair_dormancy_bonus",  "Synergy: repair×dormancy",  [0.0, 0.05, 0.15, 0.25, 0.50, 0.75, 1.0]),
    SweepParam("synergy_stress_tolerance_bonus", "Synergy: stress×tolerance", [0.0, 0.05, 0.10, 0.20, 0.40, 0.65, 1.0]),

    # --- Population / Clearance Thresholds ---
    SweepParam("carrying_capacity",   "Carrying capacity",          [1e6, 1e7, 5e7, 5e8, 2e9, 1e10, 1e11], log_scale=True),
    SweepParam("min_population",      "Min population (extinction)", [1.0, 10.0, 50.0, 100.0, 500.0, 2000.0, 10000.0], log_scale=True),
    SweepParam("clearance_threshold", "Clearance threshold",        [50.0, 200.0, 500.0, 1000.0, 5000.0, 20000.0, 100000.0], log_scale=True),
]
# fmt: on

_SWEEP_BY_FIELD: dict[str, SweepParam] = {p.field: p for p in _SWEEP_PARAMS}

_SCENARIOS: list[tuple[str, str, str, float]] = [
    # (scenario_id, abx_class, dose_level, adherence)
    ("no_abx", "none", "std", 0.70),
    ("beta_lactam", "beta_lactam", "std", 0.70),
]

_IMMUNE_STRENGTH = 0.75
_RESISTANT_FRACTION_INIT = 0.90
_DOMINANT_GENOTYPE = "R2"
_INITIAL_POPULATION = 1e6


# ---------------------------------------------------------------------------
# Core simulation helpers
# ---------------------------------------------------------------------------


def _run_episode(
    config: SimulationConfig,
    *,
    abx_class: str,
    dose_level: str,
    adherence: float,
    immune_strength: float,
    n_days: int,
    seed: int,
) -> pd.DataFrame:
    """Run one carrier episode for n_days; return a DataFrame of daily metrics."""
    rng = np.random.default_rng(seed)
    population = StrainPopulation.create_initial(
        resistant_fraction=_RESISTANT_FRACTION_INIT,
        dominant_genotype=_DOMINANT_GENOTYPE,
        initial_population=_INITIAL_POPULATION,
        rng=rng,
        strain_namespace=f"sens_{seed}",
    )

    rows: list[dict[str, Any]] = []
    for day in range(1, n_days + 1):
        population, _ = simulate_day(
            population=population,
            abx_class=abx_class,
            dose_level=dose_level,
            adherence=adherence,
            immune_strength=immune_strength,
            config=config,
            seed=seed * 100_000 + day,
        )
        total_pop = float(population.total_population)
        rf = float(compute_resistant_fraction(population.genomes, population.populations))
        p_clear = float(compute_clearance_probability(population, immune_strength, config))
        n_strains = int(population.genomes.shape[0])
        rows.append(
            {
                "day": day,
                "seed": seed,
                "resistant_fraction": rf,
                "total_population": total_pop,
                "p_clearance": p_clear,
                "n_strains": n_strains,
            }
        )
        if total_pop < config.min_population:
            # Episode cleared – fill remaining days with cleared state
            for d in range(day + 1, n_days + 1):
                rows.append(
                    {
                        "day": d,
                        "seed": seed,
                        "resistant_fraction": 0.0,
                        "total_population": 0.0,
                        "p_clearance": 1.0,
                        "n_strains": 0,
                    }
                )
            break

    return pd.DataFrame(rows)


def _run_param_sweep(
    base_config: SimulationConfig,
    sweep: SweepParam,
    scenario_id: str,
    abx_class: str,
    dose_level: str,
    adherence: float,
    n_days: int,
    n_seeds: int,
    seed_offset: int,
) -> pd.DataFrame:
    """Sweep one parameter across its defined values; return all daily rows."""
    all_rows: list[pd.DataFrame] = []
    for value in sweep.values:
        config = replace(base_config, **{sweep.field: value})
        for seed in range(n_seeds):
            df = _run_episode(
                config,
                abx_class=abx_class,
                dose_level=dose_level,
                adherence=adherence,
                immune_strength=_IMMUNE_STRENGTH,
                n_days=n_days,
                seed=seed_offset + seed,
            )
            df["param"] = sweep.field
            df["param_value"] = value
            df["scenario"] = scenario_id
            all_rows.append(df)
    return pd.concat(all_rows, ignore_index=True)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

_METRICS = [
    ("resistant_fraction", "Resistant fraction", (0, 1)),
    ("total_population", "Total population", None),
    ("p_clearance", "p_clearance / day", (0, 1)),
    ("n_strains", "Active strains", None),
]


def _value_colormap(values: list[float]) -> list[Any]:
    cmap = plt.cm.plasma
    n = len(values)
    return [cmap(i / max(n - 1, 1)) for i in range(n)]


def _plot_param_figure(
    sweep: SweepParam,
    traj: pd.DataFrame,
    scenarios: list[str],
    figures_dir: Path,
) -> None:
    n_scenarios = len(scenarios)
    n_metrics = len(_METRICS)
    fig, axes = plt.subplots(
        n_metrics,
        n_scenarios,
        figsize=(6 * n_scenarios, 4 * n_metrics),
        squeeze=False,
    )
    fig.suptitle(f"Sensitivity: {sweep.label}", fontsize=14, fontweight="bold")

    colors = _value_colormap(sweep.values)

    for col, scenario_id in enumerate(scenarios):
        sub = traj[traj["scenario"] == scenario_id]
        for row, (metric, metric_label, ylim) in enumerate(_METRICS):
            ax = axes[row][col]
            for i, value in enumerate(sweep.values):
                vs = sub[sub["param_value"] == value]
                if vs.empty:
                    continue
                mean = vs.groupby("day")[metric].mean()
                std = vs.groupby("day")[metric].std().fillna(0)
                days = mean.index.values
                label = f"{value:.3g}"
                ax.plot(days, mean.values, color=colors[i], label=label, linewidth=1.8)
                ax.fill_between(
                    days,
                    (mean - std).values,
                    (mean + std).values,
                    color=colors[i],
                    alpha=0.15,
                )
            ax.set_xlabel("Day")
            ax.set_ylabel(metric_label)
            scenario_title = "no ABX" if scenario_id == "no_abx" else scenario_id.replace("_", " ")
            ax.set_title(f"{scenario_title}")
            if ylim:
                ax.set_ylim(*ylim)
            if metric == "total_population":
                ax.set_yscale("log")
                ax.set_ylim(bottom=1.0)
            if row == 0 and col == n_scenarios - 1:
                ax.legend(
                    title=sweep.label,
                    loc="upper right",
                    fontsize=7,
                    title_fontsize=7,
                )

    fig.tight_layout()
    out = figures_dir / f"sensitivity_{sweep.field}.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.name}")


def _plot_overview_figure(
    summary: pd.DataFrame,
    sweep_params: list[SweepParam],
    figures_dir: Path,
) -> None:
    """Bar chart: range of final-day resistant_fraction across sweep values, per scenario."""
    scenario_ids = summary["scenario"].unique().tolist()
    n_scenarios = len(scenario_ids)

    fig, axes = plt.subplots(
        1, n_scenarios, figsize=(max(10, len(sweep_params) * 0.7), 6), squeeze=False
    )
    fig.suptitle("Sensitivity overview – final-day resistant_fraction range", fontsize=13)

    param_labels = [p.label for p in sweep_params]
    param_fields = [p.field for p in sweep_params]

    for col, scenario_id in enumerate(scenario_ids):
        ax = axes[0][col]
        sub = summary[
            (summary["scenario"] == scenario_id) & (summary["metric"] == "resistant_fraction")
        ]
        ranges = []
        for field in param_fields:
            ps = sub[sub["param"] == field]
            if ps.empty:
                ranges.append(0.0)
                continue
            ranges.append(float(ps["mean"].max() - ps["mean"].min()))
        x = np.arange(len(param_labels))
        ax.bar(x, ranges, color="steelblue", edgecolor="white", linewidth=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(param_labels, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Range of mean resistant_fraction")
        scenario_title = "no ABX" if scenario_id == "no_abx" else scenario_id.replace("_", " ")
        ax.set_title(f"Scenario: {scenario_title}")
        ax.set_ylim(0, 1)
        ax.axhline(0.05, color="gray", linewidth=0.8, linestyle="--", label="5% threshold")
        ax.legend(fontsize=8)

    fig.tight_layout()
    out = figures_dir / "sensitivity_overview.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.name}")


# ---------------------------------------------------------------------------
# CLI entrypoint
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Micro-parameter sensitivity analysis (Calibration 7)."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=PROJECT_ROOT / "config" / "cal_micro_sensitivity.yml",
        help="Base config YAML. Only the micro: block is used.",
    )
    parser.add_argument(
        "--params",
        nargs="+",
        default=None,
        metavar="FIELD",
        help=(
            "Subset of parameters to analyse. "
            "If omitted, all calibratable parameters are swept. "
            f"Available: {', '.join(p.field for p in _SWEEP_PARAMS)}"
        ),
    )
    parser.add_argument(
        "--n-days",
        type=int,
        default=255,
        help="Days per episode (default 255 ≈ median MRSA carriage duration).",
    )
    parser.add_argument(
        "--n-seeds", type=int, default=16, help="Seeds per (param, value, scenario)."
    )
    parser.add_argument(
        "--no-abx-only",
        action="store_true",
        help="Run only the no-ABX scenario (faster).",
    )
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    # --- resolve parameter list ---
    if args.params:
        unknown = [p for p in args.params if p not in _SWEEP_BY_FIELD]
        if unknown:
            parser.error(f"Unknown parameter(s): {', '.join(unknown)}")
        sweep_params = [_SWEEP_BY_FIELD[f] for f in args.params]
    else:
        sweep_params = list(_SWEEP_PARAMS)

    scenarios = [_SCENARIOS[0]] if args.no_abx_only else _SCENARIOS

    # --- load base config ---
    settings = load_coupled_settings(args.config)
    base_config: SimulationConfig = settings.micro

    # --- output directories ---
    output_dir = args.output_dir or (
        PROJECT_ROOT / "outputs" / (datetime.now().strftime("%Y%m%d_%H%M%S") + "_MicroSensitivity")
    )
    data_dir = output_dir / "data"
    figures_dir = output_dir / "figures"
    data_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    (
        len(sweep_params)
        * sum(len(p.values) for p in sweep_params)
        // len(sweep_params)
        * len(scenarios)
        * args.n_seeds
    )
    print("=" * 72)
    print("Micro-parameter sensitivity analysis (Calibration 7)")
    print(f"  Config:     {args.config}")
    print(f"  Parameters: {len(sweep_params)}")
    print(f"  Scenarios:  {[s[0] for s in scenarios]}")
    print(f"  Days/run:   {args.n_days}   Seeds/run: {args.n_seeds}")
    est_runs = len(sweep_params) * 5 * len(scenarios) * args.n_seeds
    print(f"  ~{est_runs} episode runs × {args.n_days} days")
    print("=" * 72)

    # --- main sweep loop ---
    all_traj: list[pd.DataFrame] = []
    for sweep in sweep_params:
        print(f"\n[{sweep.field}]")
        for scenario_id, abx_class, dose_level, adherence in scenarios:
            print(f"  scenario={scenario_id} ...", end=" ", flush=True)
            df = _run_param_sweep(
                base_config,
                sweep,
                scenario_id=scenario_id,
                abx_class=abx_class,
                dose_level=dose_level,
                adherence=adherence,
                n_days=args.n_days,
                n_seeds=args.n_seeds,
                seed_offset=0,
            )
            all_traj.append(df)
            print(f"done ({len(df)} rows)")

    traj_df = pd.concat(all_traj, ignore_index=True)

    # --- summary: final-day means per (param, param_value, scenario) ---
    final_day_df = traj_df[traj_df["day"] == args.n_days].copy()
    summary_rows: list[dict[str, Any]] = []
    for metric, _, _ in _METRICS:
        agg = (
            final_day_df.groupby(["param", "param_value", "scenario"])[metric]
            .agg(["mean", "std"])
            .reset_index()
        )
        agg["metric"] = metric
        summary_rows.append(agg)
    summary_df = pd.concat(summary_rows, ignore_index=True)

    # --- save data ---
    traj_path = data_dir / "sensitivity_trajectories.parquet"
    summary_path = data_dir / "sensitivity_summary.parquet"
    traj_df.to_parquet(traj_path, index=False)
    summary_df.to_parquet(summary_path, index=False)
    print(f"\nSaved trajectories -> {traj_path}")
    print(f"Saved summary      -> {summary_path}")

    # --- per-parameter trajectory plots ---
    print("\nGenerating figures...")
    scenario_ids = [s[0] for s in scenarios]
    for sweep in sweep_params:
        sub = traj_df[traj_df["param"] == sweep.field]
        _plot_param_figure(sweep, sub, scenario_ids, figures_dir)

    # --- overview bar chart ---
    _plot_overview_figure(summary_df, sweep_params, figures_dir)

    print("=" * 72)
    print(f"All outputs in: {output_dir}")
    print("=" * 72)


if __name__ == "__main__":
    main()
