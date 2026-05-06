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
    lambda_target_lo, lambda_target_hi = 0.0046, 0.0054

    days_per_run = combined["day"].max()
    lambda_0_theoretical = beta_eff * n_car_0 / n_total
    agg = (
        combined.groupby("day")[["carriers", "lambda_t"]]
        .agg(
            carriers_mean=("carriers", "mean"),
            carriers_p5=("carriers", lambda x: x.quantile(0.05)),
            carriers_p95=("carriers", lambda x: x.quantile(0.95)),
            lambda_mean=("lambda_t", "mean"),
            lambda_p5=("lambda_t", lambda x: x.quantile(0.05)),
            lambda_p95=("lambda_t", lambda x: x.quantile(0.95)),
        )
        .reset_index()
    )

    day1_cases = combined[combined["day"] == 1]["new_cases"].values

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
    sample_ids = np.random.choice(n_runs, size=min(100, n_runs), replace=False)
    for run_id in sample_ids:
        g = combined[combined["run_id"] == run_id]
        ax.plot(g["day"], g["carriers"] / n_total * 100, color="firebrick", alpha=0.06, lw=0.6)
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

    # --- [0,1] Kalibrierung: λ(0) als Funktion von β₀ ---
    ax = axes[0, 1]
    beta_range = np.linspace(max(0.01, beta_0 * 0.3), beta_0 * 2.0, 300)
    lambda_range = beta_range * contact_attempts * (1.0 - hygiene) * n_car_0 / n_total
    ax.plot(beta_range, lambda_range, color="steelblue", lw=2, label="λ(0) = β₀ · c · (1−H) · I₀/N")
    ax.axhspan(
        lambda_target_lo,
        lambda_target_hi,
        alpha=0.2,
        color="green",
        label=f"Zielbereich [{lambda_target_lo}–{lambda_target_hi}] Tag⁻¹",
    )
    ax.axvline(beta_0, color="darkorange", lw=1.5, ls="--", alpha=0.6)
    ax.scatter(
        [beta_0],
        [lambda_0_theoretical],
        color="darkorange",
        s=150,
        zorder=5,
        label=f"β₀ = {beta_0:.3f} → λ(0) = {lambda_0_theoretical:.4f} Tag⁻¹",
    )
    ax.set_xlabel("Transmissionsrate β₀")
    ax.set_ylabel("Initiale Infektionskraft λ(0) [Tag⁻¹]")
    ax.set_title(f"Kalibrierung: λ(0) als Funktion von β₀ ({n_runs} Seeds bestätigt)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # --- [1,0] Verteilung Tage bis 50% Durchseuchung ---
    ax = axes[1, 0]
    days_50 = summaries["days_to_50pct"].dropna()
    n_reached = len(days_50)
    if n_reached > 0:
        bins = max(15, int(days_per_run / 2))
        ax.hist(
            days_50,
            bins=bins,
            range=(1, days_per_run),
            color="steelblue",
            alpha=0.75,
            label=f"Tage bis 50% (n={n_reached}/{n_runs})",
        )
        median_50 = days_50.median()
        p5, p95 = days_50.quantile(0.05), days_50.quantile(0.95)
        ax.axvline(median_50, color="navy", lw=2, ls="--", label=f"Median = {median_50:.0f} Tage")
        ax.axvspan(p5, p95, alpha=0.1, color="navy", label=f"P5–P95 = [{p5:.0f}–{p95:.0f}] Tage")
    else:
        ax.text(0.5, 0.5, "50% nicht erreicht", transform=ax.transAxes, ha="center", va="center")
    ax.set_title("Tage bis 50% Durchseuchung (geschlossenes System)")
    ax.set_xlabel("Tage [Tage]")
    ax.set_ylabel("Anzahl Läufe")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, axis="y")

    # --- [1,1] Erstinfektionen Tag 1 — Stochastische Validierung von λ(0) ---
    ax = axes[1, 1]
    expected_cases = lambda_0_theoretical * n_sus_0
    vals, counts = np.unique(day1_cases, return_counts=True)
    ax.bar(vals, counts / n_runs * 100, color="darkorange", alpha=0.75, label="Simuliert")
    ax.axvline(
        expected_cases,
        color="green",
        lw=2,
        ls="--",
        label=f"Erwartungswert = {expected_cases:.2f}\n(λ(0) × S₀ = {lambda_0_theoretical:.4f} × {n_sus_0})",
    )
    ax.set_title("Erstinfektionen Tag 1 — Stochastische Validierung λ(0)")
    ax.set_xlabel("Neue Fälle an Tag 1")
    ax.set_ylabel("Anteil Läufe [%]")
    ax.set_xticks(range(int(day1_cases.max()) + 2))
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, axis="y")

    plt.tight_layout()
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


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
        PROJECT_ROOT / "outputs" / (datetime.now().strftime("%Y%m%d_%H%M%S") + "_Single_Ward_Batch")
    )
    data_dir = output_dir / "data"
    plot_dir = output_dir / "plots"

    print("=" * 65)
    print(f"Batch-Kalibrierung: {n_runs} Laeufe")
    print(f"  Config:  {args.config}")
    print(f"  β₀={beta_0:.3f}  c={contact_attempts:.1f} Kontakte/Tag  H={hygiene:.2f}")
    lambda_target_lo, lambda_target_hi = 0.0046, 0.0054
    in_range = lambda_target_lo <= lambda_0 <= lambda_target_hi
    status = "OK  -- innerhalb Zielbereich" if in_range else "WARNUNG -- ausserhalb Zielbereich"
    print(
        f"  → β_eff = {beta_eff:.4f} Tag⁻¹  |  λ(0) = {lambda_0:.4f} Tag⁻¹  |  Ziel: {lambda_target_lo}–{lambda_target_hi} Tag⁻¹  |  {status}"
    )
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

        if (i + 1) % max(1, n_runs // 4) == 0 or i == 0:
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
