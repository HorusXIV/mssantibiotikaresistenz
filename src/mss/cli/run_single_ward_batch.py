"""Batch-Kalibrierung: N Laeufe der Einzelgitter-Simulation mit verschiedenen Seeds.

Fuehrt run_once() N-mal mit Seeds 0..N-1 aus, aggregiert die Ergebnisse und erstellt
zusammenfassende Visualisierungen. Keine Einzelausgaben pro Lauf.

Aufruf:
    mss-calibrate-batch
    mss-calibrate-batch --n-runs 1000
    mss-calibrate-batch --config config/simulation_single_ward.yml --n-runs 500
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mss.cli.run_single_ward_calibration import (
    DEFAULT_CONFIG,
    PROJECT_ROOT,
    _derive_transmission_params,
    _load_config,
    run_once,
)

# ---------------------------------------------------------------------------
# Batch-Visualisierung
# ---------------------------------------------------------------------------


def _plot_batch(
    combined: pd.DataFrame,
    summaries: pd.DataFrame,
    beta_0: float,
    contact_attempts: float,
    hygiene: float,
    beta_eff: float,
    n_total: int,
    n_sus_0: int,
    n_car_0: int,
    n_runs: int,
    plot_path: Path,
) -> None:
    """Erstellt 4 aggregierte Visualisierungen ueber alle Batch-Laeufe."""

    days_per_run = combined["day"].max()
    lambda_0_theoretical = beta_eff * n_car_0 / n_total

    # --- Aggregation pro Tag ---
    agg = (
        combined.groupby("day")[["carriers", "lambda_t", "cumulative_attacked_pct"]]
        .agg(
            carriers_mean=("carriers", "mean"),
            carriers_p5=("carriers", lambda x: x.quantile(0.05)),
            carriers_p95=("carriers", lambda x: x.quantile(0.95)),
            lambda_mean=("lambda_t", "mean"),
            lambda_p5=("lambda_t", lambda x: x.quantile(0.05)),
            lambda_p95=("lambda_t", lambda x: x.quantile(0.95)),
            attack_mean=("cumulative_attacked_pct", "mean"),
            attack_p5=("cumulative_attacked_pct", lambda x: x.quantile(0.05)),
            attack_p95=("cumulative_attacked_pct", lambda x: x.quantile(0.95)),
        )
        .reset_index()
    )

    # Laufender Mittelwert der finalen Angriffsrate (Konvergenzplot)
    final_attacks = summaries["final_attack_pct"].values
    running_mean = np.cumsum(final_attacks) / np.arange(1, n_runs + 1)

    fig, axes = plt.subplots(2, 2, figsize=(14, 9))

    subtitle = (
        f"β₀={beta_0:.3f} · c={contact_attempts:.1f} Kontakte/Tag · H={hygiene:.2f}"
        f" → β_eff={beta_eff:.3f} Tag⁻¹ | λ(0)={lambda_0_theoretical:.4f} Tag⁻¹"
    )
    fig.suptitle(
        f"Batch-Kalibrierung  |  {n_runs} Laeufe  |"
        f"  N={n_total} (S₀={n_sus_0}, I₀={n_car_0})\n{subtitle}",
        fontsize=11,
        fontweight="bold",
    )

    # --- [0,0] Stochastische I(t)-Kurven + Mittelwert ---
    ax = axes[0, 0]
    # Zufaellige Auswahl von 100 Laeufen als transparente Hintergrundlinien
    sample_ids = np.random.choice(n_runs, size=min(100, n_runs), replace=False)
    for run_id in sample_ids:
        g = combined[combined["run_id"] == run_id]
        ax.plot(
            g["day"],
            g["carriers"] / n_total * 100,
            color="firebrick",
            alpha=0.06,
            lw=0.6,
        )
    ax.plot(
        agg["day"],
        agg["carriers_mean"] / n_total * 100,
        color="firebrick",
        lw=2.5,
        label="Mittelwert I(t)/N",
    )
    ax.fill_between(
        agg["day"],
        agg["carriers_p5"] / n_total * 100,
        agg["carriers_p95"] / n_total * 100,
        alpha=0.25,
        color="firebrick",
        label="5.–95. Perzentil",
    )
    ax.set_title(f"Prävalenz I(t)/N — {min(100, n_runs)} Läufe + Mittelwert")
    ax.set_xlabel("Tag [Tage]")
    ax.set_ylabel("Prävalenz I(t)/N [%]")
    ax.set_ylim(0, 105)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # --- [0,1] Verteilung Ausbruchszeiten ---
    ax = axes[0, 1]
    days_50 = summaries["days_to_50pct"].dropna()
    days_100 = summaries["days_to_100pct"].dropna()
    n_reached_50 = len(days_50)
    n_reached_100 = len(days_100)

    bins = max(20, int(days_per_run / 20))
    bin_range = (1, days_per_run)

    if n_reached_50 > 0:
        ax.hist(
            days_50,
            bins=bins,
            range=bin_range,
            color="steelblue",
            alpha=0.65,
            label=f"Tage bis 50% (n={n_reached_50}/{n_runs})",
        )
    if n_reached_100 > 0:
        ax.hist(
            days_100,
            bins=bins,
            range=bin_range,
            color="firebrick",
            alpha=0.55,
            label=f"Tage bis 100% (n={n_reached_100}/{n_runs})",
        )

    for col, color, label in [
        ("days_to_50pct", "steelblue", "Median 50%"),
        ("days_to_100pct", "firebrick", "Median 100%"),
    ]:
        vals = summaries[col].dropna()
        if len(vals) > 0:
            ax.axvline(
                vals.median(),
                color=color,
                lw=2,
                ls="--",
                alpha=0.9,
                label=f"{label}={vals.median():.0f} Tage",
            )

    ax.set_title("Verteilung: Tage bis 50% / 100% infiziert")
    ax.set_xlabel("Tage [Tage]")
    ax.set_ylabel("Anzahl Läufe")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    # --- [1,0] Mittlere Force of Infection ---
    ax = axes[1, 0]
    ax.plot(
        agg["day"],
        agg["lambda_mean"],
        color="darkorange",
        lw=2.5,
        label="Mittleres λ(t)",
    )
    ax.fill_between(
        agg["day"],
        agg["lambda_p5"],
        agg["lambda_p95"],
        alpha=0.25,
        color="darkorange",
        label="5.–95. Perzentil",
    )
    ax.axhline(
        lambda_0_theoretical,
        color="darkorange",
        lw=1,
        ls="--",
        alpha=0.5,
        label=f"λ(0) theoretisch = {lambda_0_theoretical:.4f} Tag⁻¹",
    )
    ax.set_title("Force of Infection λ(t) — Mittelwert ± Konfidenzband")
    ax.set_xlabel("Tag [Tage]")
    ax.set_ylabel("Infektionshazard [Tag⁻¹]")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # --- [1,1] Monte-Carlo-Konvergenz ---
    ax = axes[1, 1]
    run_indices = np.arange(1, n_runs + 1)
    ax.plot(run_indices, running_mean, color="purple", lw=1.5, label="Laufender Mittelwert")
    final_val = running_mean[-1]
    ax.axhline(
        final_val,
        color="purple",
        lw=1,
        ls="--",
        alpha=0.6,
        label=f"Endwert = {final_val:.1f}%",
    )
    ax.set_title("Konvergenz: Mittlere finale Angriffsrate")
    ax.set_xlabel("Anzahl Läufe")
    ax.set_ylabel("Laufender Mittelwert [%]")
    ax.set_ylim(0, 105)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot gespeichert → {plot_path}")


# ---------------------------------------------------------------------------
# Batch-Einstiegspunkt
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Batch-Kalibrierung: N Laeufe der Einzelgitter-Simulation"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="Pfad zur YAML-Konfigurationsdatei",
    )
    parser.add_argument(
        "--n-runs",
        type=int,
        default=1000,
        help="Anzahl Simulationslaeufe (default: 1000)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Ausgabeverzeichnis (default: outputs/YYYYMMDD_HHMMSS)",
    )
    args = parser.parse_args()

    raw = _load_config(args.config)
    n_runs: int = args.n_runs

    run_cfg = raw.get("run", {})
    days: int = int(run_cfg.get("days", 100))
    pop_cfg = raw.get("population", {})
    n_sus_0: int = int(pop_cfg.get("susceptible_count", 8))
    n_car_0: int = int(pop_cfg.get("carrier_count", 2))
    n_total: int = n_sus_0 + n_car_0

    beta_0, contact_attempts, hygiene, beta_eff, _ = _derive_transmission_params(raw)
    lambda_0 = beta_eff * n_car_0 / n_total

    output_dir = args.output_dir or (
        PROJECT_ROOT / "outputs" / datetime.now().strftime("%Y%m%d_%H%M%S")
    )
    data_dir = output_dir / "data"
    plot_dir = output_dir / "plots"

    print("=" * 65)
    print(f"Batch-Kalibrierung: {n_runs} Laeufe")
    print(f"  Config:  {args.config}")
    print(f"  β₀={beta_0:.3f}  c={contact_attempts:.1f} Kontakte/Tag  H={hygiene:.2f}")
    print(f"  → β_eff = {beta_eff:.4f} Tag⁻¹  |  λ(0) = {lambda_0:.4f} Tag⁻¹")
    print(f"  Personen: N={n_total} (S₀={n_sus_0}, I₀={n_car_0})  |  Tage/Lauf: {days}")
    print(f"  Seeds: 0 bis {n_runs - 1} (deterministisch, reproduzierbar)")
    print("=" * 65)

    all_dfs: List[pd.DataFrame] = []
    all_metas: List[Dict[str, Any]] = []

    for i in range(n_runs):
        df, meta = run_once(raw, seed=i)
        df = df.assign(run_id=i)
        all_dfs.append(df)
        all_metas.append(meta)

        if (i + 1) % 100 == 0 or i == 0:
            pct_done = (i + 1) / n_runs * 100
            print(f"  [{pct_done:5.1f}%] Lauf {i + 1:>{len(str(n_runs))}}/{n_runs} abgeschlossen")

    print("=" * 65)
    print("Aggregiere Ergebnisse...")

    combined = pd.concat(all_dfs, ignore_index=True)
    summaries = pd.DataFrame(all_metas)

    # Zusammenfassung
    for milestone in [25, 50, 75, 100]:
        col = f"days_to_{milestone}pct"
        vals = summaries[col].dropna()
        reached = len(vals)
        if reached > 0:
            print(
                f"  {milestone:3d}% erreicht in {reached}/{n_runs} Laeufen"
                f"  |  Median: {vals.median():.0f} Tage"
                f"  |  [P5={vals.quantile(0.05):.0f}, P95={vals.quantile(0.95):.0f}]"
            )
        else:
            print(f"  {milestone:3d}% nicht erreicht in {days} Tagen (0/{n_runs} Laeufe)")

    mean_final = summaries["final_attack_pct"].mean()
    print(f"\n  Mittlere finale Angriffsrate: {mean_final:.1f}%")

    # Parquet speichern
    data_dir.mkdir(parents=True, exist_ok=True)
    combined_path = data_dir / "batch_all_runs.parquet"
    summary_path = data_dir / "batch_summary.parquet"
    combined.to_parquet(combined_path, index=False)
    summaries.to_parquet(summary_path, index=False)

    print("\nGespeichert:")
    print(f"  Alle Laeufe → {combined_path}")
    print(f"  Zusammenfassung → {summary_path}")

    # Visualisierungen
    plot_path = plot_dir / "01_single_ward_batch.png"
    _plot_batch(
        combined=combined,
        summaries=summaries,
        beta_0=beta_0,
        contact_attempts=contact_attempts,
        hygiene=hygiene,
        beta_eff=beta_eff,
        n_total=n_total,
        n_sus_0=n_sus_0,
        n_car_0=n_car_0,
        n_runs=n_runs,
        plot_path=plot_path,
    )
    print(f"  Plot  → {plot_path}")


if __name__ == "__main__":
    main()
