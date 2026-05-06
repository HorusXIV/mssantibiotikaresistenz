"""Parameter-Sweep Kalibrierung fuer die realistische Simulation.

Fuehrt eine strukturierte Sensitivitaetsanalyse durch: Ein Parameter wird ueber
ein definiertes Raster variiert, die Simulation laeuft fuer jeden Wert, und der
Effekt auf eine Zielgroesse wird geplottet. Der Wert mit dem geringsten Abstand
zum Literatur-Zielwert wird als optimaler Wert markiert.

Methodisches Prinzip: Analogon zur analytischen beta0-Kalibrierung, aber mit
der vollstaendigen stochastischen Simulation statt einer geschlossenen Formel.
Das Vorgehen entspricht einem sequentiellen Koordinatenabstieg: Parameter werden
einzeln kalibriert, alle anderen bleiben auf ihren aktuellen Schaetzwerten.

Aufruf:
    mss-sweep --sweep config/calibration/phase1_isolation_effectiveness.yml
    mss-sweep --sweep config/calibration/phase1_isolation_effectiveness.yml --n-seeds 5
    mss-sweep --sweep config/calibration/phase1_isolation_effectiveness.yml --config config/simulation_realistic.yml
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

DEFAULT_SWEEP_CONFIG = (
    PROJECT_ROOT / "config" / "calibration" / "phase1_isolation_effectiveness.yml"
)


# ---------------------------------------------------------------------------
# Datenstruktur fuer Sweep-Parameter
# ---------------------------------------------------------------------------


@dataclass
class SweepParameter:
    section: str  # z.B. "macro" oder "population.susceptible_template"
    name: str  # z.B. "base_isolation_effectiveness"
    values: list[float]  # explizite Werte ODER via min/max/n_steps erzeugt


# ---------------------------------------------------------------------------
# Hilfsfunktionen
# ---------------------------------------------------------------------------


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
    return SweepParameter(section=section, name=name, values=values)


def _set_param(raw: dict[str, Any], section: str, name: str, value: float) -> None:
    """Setzt raw[section][name] = value.

    Unterstuetzt einfache Abschnitte ("macro") und verschachtelte Pfade
    ("population.susceptible_template").
    """
    parts = section.split(".")
    node = raw
    for part in parts:
        node = node[part]
    node[name] = value


# ---------------------------------------------------------------------------
# Sweep-Logik
# ---------------------------------------------------------------------------


def run_sweep(
    raw_base: dict[str, Any],
    spec: SweepParameter,
    target_metric: str,
    target_value: float,
    n_seeds: int,
    use_micro: bool = True,
) -> pd.DataFrame:
    """Fuehrt den Sweep durch und gibt einen DataFrame mit Ergebnissen zurueck.

    Fuer jeden Rasterpunkt werden n_seeds Simulationen gemittelt. Neben dem
    Mittelwert wird auch die Standardabweichung gespeichert (Suffix '_std'),
    um die Streuung ueber Seeds im Plot darstellen zu koennen.

    Returns:
        DataFrame mit Spalten: param_value, <metrik>, <metrik>_std, ...
        fuer alle Metriken aus dem Meta-Dict.
    """
    rows = []
    for val in spec.values:
        raw = copy.deepcopy(raw_base)
        _set_param(raw, spec.section, spec.name, val)

        seed_metas = [
            run_realistic_once(raw, seed=s, use_micro=use_micro)[1] for s in range(n_seeds)
        ]
        metric_keys = [k for k in seed_metas[0] if k != "seed"]

        mean_metrics = {k: float(np.mean([m[k] for m in seed_metas])) for k in metric_keys}
        std_metrics = {
            f"{k}_std": float(np.std([m[k] for m in seed_metas], ddof=0)) for k in metric_keys
        }

        rows.append({"param_value": val, **mean_metrics, **std_metrics})
        observed = mean_metrics[target_metric]
        abstand = abs(observed - target_value)
        print(
            f"  {spec.section}.{spec.name}={val:.4f}"
            f"  →  {target_metric}={observed:.4f}"
            f"  (Abstand={abstand:.4f})"
        )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Visualisierung
# ---------------------------------------------------------------------------


def _plot_sweep(
    df: pd.DataFrame,
    spec: SweepParameter,
    target_metric: str,
    target_value: float,
    target_label: str,
    plot_path: Path,
    n_seeds: int,
) -> None:
    """Plottet Sweep-Ergebnis: Parameter-Wert vs. beobachtete Metrik."""
    best_idx = int((df[target_metric] - target_value).abs().argmin())
    best_val = float(df.loc[best_idx, "param_value"])
    best_obs = float(df.loc[best_idx, target_metric])

    std_col = f"{target_metric}_std"
    has_std = std_col in df.columns and n_seeds > 1

    fig, ax = plt.subplots(figsize=(9, 5.5))

    # Std-Band
    if has_std:
        ax.fill_between(
            df["param_value"],
            df[target_metric] - df[std_col],
            df[target_metric] + df[std_col],
            alpha=0.18,
            color="steelblue",
            label=f"±1 Std. ({n_seeds} Seeds)",
        )

    # Hauptlinie
    ax.plot(
        df["param_value"],
        df[target_metric],
        color="steelblue",
        lw=2.5,
        marker="o",
        markersize=7,
        label=target_metric,
    )

    # Zielwert-Linie
    ax.axhline(target_value, color="green", lw=1.5, ls="--", label=f"Ziel = {target_value:.3g}")

    # Optimum-Punkt + vertikale Hilfslinie
    ax.scatter([best_val], [best_obs], color="darkorange", s=180, zorder=5)
    ax.axvline(best_val, color="darkorange", lw=1.0, ls=":", alpha=0.6)

    # Optimum direkt am Punkt annotieren
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

    ax.set_xlabel(f"{spec.section}.{spec.name}", fontsize=11)
    ax.set_ylabel(target_metric, fontsize=11)
    ax.set_title(
        f"Sweep: {spec.name}  →  {target_metric}\n" f"Ziel: {target_value:.3g}  |  {target_label}",
        fontsize=10,
    )
    ax.legend(fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3)
    ax.grid(alpha=0.3)
    fig.subplots_adjust(bottom=0.20)
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI-Einstiegspunkt
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parameter-Sweep Kalibrierung fuer die realistische Simulation"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG_PATH,
        help="Pfad zur Basis-YAML-Konfiguration (default: simulation_realistic.yml)",
    )
    parser.add_argument(
        "--sweep",
        type=Path,
        default=DEFAULT_SWEEP_CONFIG,
        help="Pfad zur Sweep-Konfigurationsdatei",
    )
    parser.add_argument(
        "--n-seeds",
        type=int,
        default=None,
        help="Anzahl Seeds zur Mittelung pro Rasterpunkt (ueberschreibt Wert in YAML)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Ausgabeverzeichnis (default: outputs/YYYYMMDD_HHMMSS_Sweep_<name>)",
    )
    args = parser.parse_args()

    sweep_cfg = _load_sweep_config(args.sweep)
    raw_base = _load_base_config(args.config)

    # Optionaler run.days-Override aus der Sweep-YAML (fuer schnellere Phase-1-Sweeps)
    if "run_days_override" in sweep_cfg:
        raw_base.setdefault("run", {})["days"] = int(sweep_cfg["run_days_override"])

    spec = _parse_sweep_parameter(sweep_cfg["parameter"])
    target_cfg = sweep_cfg["target"]
    target_metric: str = str(target_cfg["metric"])
    target_value: float = float(target_cfg["value"])
    target_label: str = str(target_cfg.get("label", "Literatur"))
    n_seeds: int = args.n_seeds if args.n_seeds is not None else int(sweep_cfg.get("n_seeds", 3))
    use_micro: bool = bool(sweep_cfg.get("use_micro", True))

    output_dir = args.output_dir or (
        PROJECT_ROOT
        / "outputs"
        / (datetime.now().strftime("%Y%m%d_%H%M%S") + f"_Sweep_{spec.name}")
    )
    data_dir = output_dir / "data"
    plot_dir = output_dir / "plots"

    run_days = raw_base.get("run", {}).get("days", 365)

    print("=" * 65)
    print(f"Parameter-Sweep: {spec.section}.{spec.name}")
    print(f"  Werte: [{', '.join(f'{v:.3f}' for v in spec.values)}]")
    print(f"  Zielgroesse: {target_metric} = {target_value}  ({target_label})")
    print(f"  Seeds pro Rasterpunkt: {n_seeds}")
    print(f"  Simulationsdauer: {run_days} Tage")
    print(f"  Mikro-Simulator: {'aktiviert' if use_micro else 'deaktiviert (use_micro: false)'}")
    print(f"  Config: {args.config}")
    if run_days > 180:
        print(
            "  HINWEIS: run.days > 180 — fuer schnellere Phase-1/2-Sweeps empfiehlt"
            " sich run_days_override: 90 in der Sweep-YAML."
        )
    print("=" * 65)

    df = run_sweep(
        raw_base=raw_base,
        spec=spec,
        target_metric=target_metric,
        target_value=target_value,
        n_seeds=n_seeds,
        use_micro=use_micro,
    )

    best_idx = int((df[target_metric] - target_value).abs().argmin())
    best_val = float(df.loc[best_idx, "param_value"])
    best_obs = float(df.loc[best_idx, target_metric])

    print("=" * 65)
    print(
        f"BESTER WERT: {spec.section}.{spec.name} = {best_val:.4f}"
        f"  →  {target_metric} = {best_obs:.4f}"
        f"  (Abstand zum Ziel: {abs(best_obs - target_value):.4f})"
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
    )

    print("\nGespeichert:")
    print(f"  Daten → {data_path}")
    print(f"  Plot  → {plot_path}")
    print("=" * 65)


if __name__ == "__main__":
    main()
