"""Parameter-sweep calibration for the realistic simulation.

Structured sensitivity analysis: one parameter is varied over a defined grid,
the simulation runs for each value, and the effect on a target metric is plotted.
The value closest to the literature target is marked as optimal.

Analogous to the analytic beta_0 calibration, but using the full stochastic
simulation instead of a closed-form formula. Follows a sequential coordinate
descent: parameters are calibrated one at a time while the rest stay at their
current estimates.

Usage:
    mss-sweep --sweep config/calibration/cal2_proximity_decay.yml
    mss-sweep --sweep config/calibration/cal2_proximity_decay.yml --n-seeds 5
    mss-sweep --sweep config/calibration/cal2_proximity_decay.yml --config config/simulation_realistic.yml
"""

from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from mss.cli.run_coupled_simulation import PROJECT_ROOT, DEFAULT_CONFIG_PATH, run_realistic_once

DEFAULT_SWEEP_CONFIG = PROJECT_ROOT / "config" / "calibration" / "cal2_proximity_decay.yml"


@dataclass
class SweepParameter:
    section: str  # e.g. "macro" or "population.susceptible_template"
    name: str  # e.g. "base_isolation_effectiveness"
    values: list[float]  # explicit values, or generated from min/max/n_steps
    unit: str = ""  # axis unit for the plot, e.g. "1/grid cell"


def _load_sweep_config(path: Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _load_base_config(path: Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _parse_sweep_parameter(param_cfg: dict[str, Any]) -> SweepParameter:
    section = str(param_cfg["section"])
    name = str(param_cfg["name"])
    if "values" in param_cfg:
        values = [float(v) for v in param_cfg["values"]]
    else:
        lo = float(param_cfg["min"])
        hi = float(param_cfg["max"])
        n = int(param_cfg.get("n_steps", 7))
        values = list(np.linspace(lo, hi, n))
    unit = str(param_cfg.get("unit", ""))
    return SweepParameter(section=section, name=name, values=values, unit=unit)


def _set_param(raw: dict[str, Any], section: str, name: str, value: float) -> None:
    """Set raw[section][name] = value.

    Supports simple sections ("macro") and nested paths
    ("population.susceptible_template").
    """
    parts = section.split(".")
    node = raw
    for part in parts:
        node = node[part]
    node[name] = value


def run_sweep(
    raw_base: dict[str, Any],
    spec: SweepParameter,
    target_metric: str,
    target_value: float,
    n_seeds: int,
    use_micro: bool = True,
    baseline_overrides: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """Run the sweep and return a DataFrame of results.

    Each grid point averages n_seeds simulations; the standard deviation is
    stored next to the mean (suffix '_std') to show the seed spread in the plot.

    If ``baseline_overrides`` is set, a per-seed counterfactual run (e.g. isolation
    off) is computed and the derived metric
    ``acquisition_reduction = 1 - n_transmissions_sweep / n_transmissions_baseline``
    (paired per seed) is added. The ratio is beta_0-independent because the level
    cancels out. Basis is the cumulative in-hospital transmission count over the
    whole run (high event count, low noise).

    Returns:
        DataFrame with columns param_value, <metric>, <metric>_std, ... for every
        metric in the meta dict.
    """
    baseline_acq: list[float] | None = None
    if baseline_overrides:
        raw_bl = copy.deepcopy(raw_base)
        for section, overrides in baseline_overrides.items():
            node = raw_bl.setdefault(section, {})
            for key, value in overrides.items():
                node[key] = value
        baseline_acq = [
            float(run_realistic_once(raw_bl, seed=s, use_micro=use_micro)[1]["n_transmissions"])
            for s in range(n_seeds)
        ]
        print(f"  Baseline (counterfactual) transmissions off = {float(np.mean(baseline_acq)):.1f}")

    rows = []
    for val in spec.values:
        raw = copy.deepcopy(raw_base)
        _set_param(raw, spec.section, spec.name, val)

        seed_metas = [
            run_realistic_once(raw, seed=s, use_micro=use_micro)[1] for s in range(n_seeds)
        ]
        if baseline_acq is not None:
            for s, meta in enumerate(seed_metas):
                a_off = baseline_acq[s]
                meta["acquisition_reduction"] = (
                    1.0 - float(meta["n_transmissions"]) / a_off if a_off > 0 else 0.0
                )
        metric_keys = [k for k in seed_metas[0] if k != "seed"]

        mean_metrics = {k: float(np.mean([m[k] for m in seed_metas])) for k in metric_keys}
        std_metrics = {
            f"{k}_std": float(np.std([m[k] for m in seed_metas], ddof=0)) for k in metric_keys
        }

        rows.append({"param_value": val, **mean_metrics, **std_metrics})
        observed = mean_metrics[target_metric]
        distance = abs(observed - target_value)
        print(
            f"  {spec.section}.{spec.name}={val:.4f}"
            f"  ->  {target_metric}={observed:.4f}"
            f"  (distance={distance:.4f})"
        )

    return pd.DataFrame(rows)


def _plot_sweep(
    df: pd.DataFrame,
    spec: SweepParameter,
    target_metric: str,
    target_value: float,
    target_label: str,
    plot_path: Path,
    n_seeds: int,
    target_unit: str = "",
) -> None:
    """Plot the sweep result: parameter value vs. observed metric."""
    x_label = (
        f"{spec.section}.{spec.name} [{spec.unit}]" if spec.unit else f"{spec.section}.{spec.name}"
    )
    y_label = f"{target_metric} [{target_unit}]" if target_unit else target_metric
    best_idx = int((df[target_metric] - target_value).abs().argmin())
    best_val = float(df.loc[best_idx, "param_value"])
    best_obs = float(df.loc[best_idx, target_metric])

    std_col = f"{target_metric}_std"
    has_std = std_col in df.columns and n_seeds > 1

    fig, ax = plt.subplots(figsize=(9, 5.5))

    # Standard-deviation band across seeds
    if has_std:
        ax.fill_between(
            df["param_value"],
            df[target_metric] - df[std_col],
            df[target_metric] + df[std_col],
            alpha=0.18,
            color="steelblue",
            label=f"±1 std ({n_seeds} seeds)",
        )

    ax.plot(
        df["param_value"],
        df[target_metric],
        color="steelblue",
        lw=2.5,
        marker="o",
        markersize=7,
        label=target_metric,
    )

    ax.axhline(target_value, color="green", lw=1.5, ls="--", label=f"target = {target_value:.3g}")

    # Optimum point with vertical guide line
    ax.scatter([best_val], [best_obs], color="darkorange", s=180, zorder=5)
    ax.axvline(best_val, color="darkorange", lw=1.0, ls=":", alpha=0.6)

    ax.annotate(
        f"{spec.name} = {best_val:.4f}\n{target_metric} = {best_obs:.4f}",
        xy=(best_val, best_obs),
        xytext=(12, 12),
        textcoords="offset points",
        fontsize=8,
        color="darkorange",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="darkorange", alpha=0.85),
        arrowprops=dict(arrowstyle="->", color="darkorange", lw=0.8),
    )

    ax.set_xlabel(x_label, fontsize=11)
    ax.set_ylabel(y_label, fontsize=11)
    ax.set_title(
        f"Sweep: {spec.name}  ->  {target_metric}\n"
        f"target: {target_value:.3g}  |  {target_label}",
        fontsize=10,
    )
    ax.legend(fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3)
    ax.grid(alpha=0.3)
    fig.subplots_adjust(bottom=0.20)
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parameter-sweep calibration for the realistic simulation"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG_PATH,
        help="Path to the base YAML config (default: simulation_realistic.yml)",
    )
    parser.add_argument(
        "--sweep",
        type=Path,
        default=DEFAULT_SWEEP_CONFIG,
        help="Path to the sweep config file",
    )
    parser.add_argument(
        "--n-seeds",
        type=int,
        default=None,
        help="Number of seeds averaged per grid point (overrides the YAML value)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory (default: outputs/YYYYMMDD_HHMMSS_Sweep_<name>)",
    )
    args = parser.parse_args()

    sweep_cfg = _load_sweep_config(args.sweep)
    raw_base = _load_base_config(args.config)

    # Optional run.days override from the sweep YAML (for faster macro sweeps)
    if "run_days_override" in sweep_cfg:
        raw_base.setdefault("run", {})["days"] = int(sweep_cfg["run_days_override"])

    # Optional base_overrides: set base-config fields for the sweep, e.g. disable
    # antibiotics to calibrate one parameter in isolation.
    for section, overrides in sweep_cfg.get("base_overrides", {}).items():
        node = raw_base.setdefault(section, {})
        for key, value in overrides.items():
            node[key] = value

    spec = _parse_sweep_parameter(sweep_cfg["parameter"])
    target_cfg = sweep_cfg["target"]
    target_metric: str = str(target_cfg["metric"])
    target_value: float = float(target_cfg["value"])
    target_label: str = str(target_cfg.get("label", "literature"))
    target_unit: str = str(target_cfg.get("unit", ""))
    n_seeds: int = args.n_seeds if args.n_seeds is not None else int(sweep_cfg.get("n_seeds", 3))
    use_micro: bool = bool(sweep_cfg.get("use_micro", True))
    baseline_overrides: dict[str, Any] = sweep_cfg.get("baseline_overrides", {})

    output_dir = args.output_dir or (
        PROJECT_ROOT
        / "outputs"
        / (datetime.now().strftime("%Y%m%d_%H%M%S") + f"_Sweep_{spec.name}")
    )
    data_dir = output_dir / "data"
    plot_dir = output_dir / "plots"

    run_days = raw_base.get("run", {}).get("days", 365)

    print("=" * 65)
    print(f"Parameter sweep: {spec.section}.{spec.name}")
    print(f"  Values: [{', '.join(f'{v:.3f}' for v in spec.values)}]")
    print(f"  Target metric: {target_metric} = {target_value}  ({target_label})")
    print(f"  Seeds per grid point: {n_seeds}")
    print(f"  Simulated days: {run_days}")
    print(f"  Micro simulator: {'on' if use_micro else 'off (use_micro: false)'}")
    print(f"  Config: {args.config}")
    if run_days > 180:
        print(
            "  NOTE: run.days > 180; for faster macro sweeps set"
            " run_days_override: 90 in the sweep YAML."
        )
    print("=" * 65)

    df = run_sweep(
        raw_base=raw_base,
        spec=spec,
        target_metric=target_metric,
        target_value=target_value,
        n_seeds=n_seeds,
        use_micro=use_micro,
        baseline_overrides=baseline_overrides,
    )

    best_idx = int((df[target_metric] - target_value).abs().argmin())
    best_val = float(df.loc[best_idx, "param_value"])
    best_obs = float(df.loc[best_idx, target_metric])

    print("=" * 65)
    print(
        f"BEST VALUE: {spec.section}.{spec.name} = {best_val:.4f}"
        f"  ->  {target_metric} = {best_obs:.4f}"
        f"  (distance to target: {abs(best_obs - target_value):.4f})"
    )

    data_dir.mkdir(parents=True, exist_ok=True)
    data_path = data_dir / f"sweep_{spec.name}.parquet"
    df.to_parquet(data_path, index=False)

    plot_path = plot_dir / f"sweep_{spec.name}.png"
    _plot_sweep(
        df=df,
        spec=spec,
        target_metric=target_metric,
        target_value=target_value,
        target_label=target_label,
        plot_path=plot_path,
        n_seeds=n_seeds,
        target_unit=target_unit,
    )

    print("\nSaved:")
    print(f"  Data -> {data_path}")
    print(f"  Plot -> {plot_path}")
    print("=" * 65)


if __name__ == "__main__":
    main()
