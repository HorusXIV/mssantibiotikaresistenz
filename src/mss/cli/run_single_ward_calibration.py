"""Einzelgitter-Kalibrierung fuer die Makro-Simulation.

Dieses Skript laeuft den MacroSimulator auf einem einzigen, geschlossenen
Gitter (eine Ward, N Personen) ohne Mikro-Simulation.

Ziel: Realitaetsnahe Werte fuer die wichtigsten Transmissionsparameter finden,
indem beobachtet wird wie lange es dauert bis alle Personen infiziert sind.

Das Modell ist ein diskretes stochastisches SI-Modell ohne Erholung:
    lambda(t) = beta_eff * I(t) / N         [Tag^-1]
    beta_eff  = base_transmission_rate * daily_contact_attempts * (1 - base_hygiene)  [Tag^-1]
    p_inf     = 1 - exp(-lambda(t))          [dimensionslos]

Aufruf (Einzellauf):
    mss-calibrate
    mss-calibrate --config config/calibration/cal1_simulation_single_ward.yml

Aufruf (Ensemble über mehrere Seeds):
    mss-calibrate --n-runs 1000
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CONFIG = PROJECT_ROOT / "config" / "calibration" / "cal1_simulation_single_ward.yml"


# ---------------------------------------------------------------------------
# Config laden
# ---------------------------------------------------------------------------


def _load_config(config_path: Path) -> Dict[str, Any]:
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _build_macro_config(raw: Dict[str, Any]):
    from mss.simulation.macro.config import SimulationConfig

    macro_raw = raw.get("macro", {})
    return SimulationConfig(**{k: v for k, v in macro_raw.items() if hasattr(SimulationConfig, k)})


def _derive_transmission_params(raw: Dict[str, Any]) -> Tuple[float, float, float, float, float]:
    """Gibt (beta_0, contact_attempts, hygiene, beta_eff, n_car_0/n_total) zurueck."""
    pop = raw.get("population", {})
    macro_raw = raw.get("macro", {})
    n_sus_0 = int(pop.get("susceptible_count", 8))
    n_car_0 = int(pop.get("carrier_count", 2))
    n_total = n_sus_0 + n_car_0
    beta_0 = float(macro_raw.get("base_transmission_rate", 0.09))
    contact_attempts = float(macro_raw.get("daily_contact_attempts", 18.0))
    hygiene = float(macro_raw.get("base_hygiene", 0.65))
    beta_eff = beta_0 * contact_attempts * (1.0 - hygiene)
    return beta_0, contact_attempts, hygiene, beta_eff, n_car_0 / n_total


# ---------------------------------------------------------------------------
# Patienten erstellen
# ---------------------------------------------------------------------------


def _create_patients(raw: Dict[str, Any]) -> List[Any]:
    from mss.domain.patient import Department, HealthState, Patient

    pop = raw.get("population", {})
    n_sus = int(pop.get("susceptible_count", 8))
    n_car = int(pop.get("carrier_count", 2))

    patients = []

    for i in range(n_sus):
        p = Patient(
            patient_id=f"sus_{i:03d}",
            state=HealthState.SUSCEPTIBLE,
            p_clearance=0.0,  # Keine spontane Selbstheilung
            department=Department.WARD,
        )
        patients.append(p)

    for i in range(n_car):
        p = Patient(
            patient_id=f"car_{i:03d}",
            state=HealthState.CARRIER,
            episode_id=f"ep_initial_{i:03d}",
            p_clearance=0.0,  # Keine spontane Selbstheilung
            dominant_genotype="S",
            resistant_fraction=0.0,
            department=Department.WARD,
        )
        patients.append(p)

    return patients


# ---------------------------------------------------------------------------
# Kernsimulation (wiederverwendbar fuer Batch-Laeufe)
# ---------------------------------------------------------------------------


def run_once(raw: Dict[str, Any], seed: int) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Fuehrt einen einzelnen Kalibrierungslauf durch und gibt (df, metadata) zurueck.

    Laeuft immer fuer die volle konfigurierte Laufzeit (kein Fruehausstieg), damit
    Batch-Laeufe zeitlich vergleichbare Zeitreihen liefern.

    Parameters
    ----------
    raw : dict
        Geladene YAML-Konfiguration.
    seed : int
        Zufallsseed fuer diesen Lauf.

    Returns
    -------
    df : pd.DataFrame
        Taeglich: day, susceptible, carriers, new_cases, prevalence_pct,
        cumulative_attacked_pct, lambda_t.
    metadata : dict
        seed, days_to_{25,50,75,100}pct (int oder None), final_attack_pct (float).
    """
    from mss.domain.patient import Department, HealthState
    from mss.simulation.macro.simulator import MacroSimulator

    run_cfg = raw.get("run", {})
    days: int = int(run_cfg.get("days", 100))
    run_id: str = str(run_cfg.get("run_id", "single_ward_calibration"))

    pop_cfg = raw.get("population", {})
    n_sus_0: int = int(pop_cfg.get("susceptible_count", 8))
    n_car_0: int = int(pop_cfg.get("carrier_count", 2))
    n_total: int = n_sus_0 + n_car_0

    _, _, _, beta_eff, _ = _derive_transmission_params(raw)

    macro_cfg = _build_macro_config(raw)
    macro = MacroSimulator(config=macro_cfg, n_hospitals=1, seed=seed)

    hospital_id = "hospital_001"
    for patient in _create_patients(raw):
        macro.admit(patient, hospital_id, Department.WARD)

    records: List[Dict] = []
    milestones: Dict[int, int] = {}
    prev_carriers = n_car_0
    cumulative = n_car_0  # initial carriers zaehlen als "je infiziert"

    for day in range(1, days + 1):
        # Zustand zu Tagesbeginn erfassen
        current_before = macro.get_patients(hospital_id)
        n_c_before = sum(1 for p in current_before if p.state == HealthState.CARRIER)
        lambda_t = beta_eff * n_c_before / n_total

        # Tagesupdate ausfuehren
        macro.step(micro_simulator=None, run_id=run_id)

        # _colonize() setzt p_clearance=0.02 bei jeder Neuinfektion (designed fuer Mikro,
        # das den Wert direkt ueberschreibt). Ohne Mikro bleibt 0.02 bestehen -> Oszillation.
        # Fix: nach jedem Step alle Patienten auf p_clearance=0.0 zuruecksetzen.
        for p in macro.get_patients(hospital_id):
            p.p_clearance = 0.0

        # Zustand nach Tagesupdate erfassen
        current_after = macro.get_patients(hospital_id)
        n_c = sum(1 for p in current_after if p.state == HealthState.CARRIER)
        new_cases = max(0, n_c - prev_carriers)
        cumulative += new_cases

        prevalence_pct = n_c / n_total * 100
        cumulative_pct = cumulative / n_total * 100

        for milestone in [25, 50, 75, 100]:
            if milestone not in milestones and cumulative_pct >= milestone:
                milestones[milestone] = day

        records.append(
            {
                "day": day,
                "susceptible": n_total - n_c,
                "carriers": n_c,
                "new_cases": new_cases,
                "prevalence_pct": prevalence_pct,
                "cumulative_attacked_pct": cumulative_pct,
                "lambda_t": lambda_t,
            }
        )
        prev_carriers = n_c

    df = pd.DataFrame(records)
    metadata: Dict[str, Any] = {
        "seed": seed,
        "days_to_25pct": milestones.get(25),
        "days_to_50pct": milestones.get(50),
        "days_to_75pct": milestones.get(75),
        "days_to_100pct": milestones.get(100),
        "final_attack_pct": df["cumulative_attacked_pct"].iloc[-1],
    }
    return df, metadata


# ---------------------------------------------------------------------------
# Einzellauf-Visualisierung
# ---------------------------------------------------------------------------


def _plot(
    df: pd.DataFrame,
    beta_eff: float,
    beta_0: float,
    contact_attempts: float,
    hygiene: float,
    n_total: int,
    n_sus_0: int,
    n_car_0: int,
    seed: int,
    plot_path: Path,
) -> None:
    lambda_target_lo, lambda_target_hi = 0.0046, 0.0054

    carriers = df["carriers"]
    new_cases = df["new_cases"]
    cumulative_pct = df["cumulative_attacked_pct"]
    days = df["day"]
    lambda_0_theoretical = beta_eff * n_car_0 / n_total
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))

    subtitle = (
        f"β₀={beta_0:.3f} · c={contact_attempts:.1f} Kontakte/Tag · H={hygiene:.2f}"
        f" → β_eff={beta_eff:.3f} Tag⁻¹ | λ(0)={lambda_0_theoretical:.4f} Tag⁻¹"
    )
    fig.suptitle(
        f"Einzelgitter-Kalibrierung  |  N={n_total} (S₀={n_sus_0}, I₀={n_car_0})  |  Seed={seed}\n{subtitle}",
        fontsize=11,
        fontweight="bold",
    )

    # --- [0,0] Kalibrierung: λ(0) als Funktion von β₀ ---
    ax = axes[0, 0]
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
    ax.set_xlabel("Transmissionsrate β₀ [pro Kontakt*Carrier]")
    ax.set_ylabel("Initiale Infektionskraft λ(0) [Tag⁻¹]")
    ax.set_title("Kalibrierung: λ(0) als Funktion von β₀")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # --- [0,1] SI-Dynamik ---
    ax = axes[0, 1]
    ax.plot(days, df["susceptible"], color="steelblue", lw=2, label="Susceptible S(t)")
    ax.plot(days, carriers, color="firebrick", lw=2, label="Carrier I(t)")
    ax.axhline(n_total, color="gray", lw=1, ls="--", alpha=0.5, label=f"Total N={n_total}")
    ax.set_title("SI-Dynamik")
    ax.set_xlabel("Tag [Tage]")
    ax.set_ylabel("Patienten [Anzahl]")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    ax.set_ylim(0, n_total * 1.15)

    # --- [1,0] Tägliche Inzidenz ---
    ax = axes[1, 0]
    ax.bar(days, new_cases, color="firebrick", alpha=0.7, label="Neue Fälle/Tag")
    ax.set_title("Tägliche Inzidenz")
    ax.set_xlabel("Tag [Tage]")
    ax.set_ylabel("Neue Fälle [Patienten]")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, axis="y")
    if n_total > 0:
        ax2 = ax.twinx()
        ax2.set_ylabel("Rate [pro 1000 Patienten-Tage]", fontsize=8, color="gray")
        ax2.set_ylim(0, max(ax.get_ylim()[1] / n_total * 1000, 0.1))
        ax2.tick_params(colors="gray")

    # --- [1,1] Kumulative Angriffsrate (geschlossenes System) ---
    ax = axes[1, 1]
    ax.plot(days, cumulative_pct, color="purple", lw=2)
    ax.fill_between(days, cumulative_pct, alpha=0.15, color="purple")
    for pct_target, label in {25: "25%", 50: "50%", 75: "75%", 100: "100%"}.items():
        hit = df[df["cumulative_attacked_pct"] >= pct_target]
        if not hit.empty:
            day_hit = hit.iloc[0]["day"]
            ax.axvline(day_hit, color="gray", lw=1, ls=":", alpha=0.7)
            ax.annotate(
                f"{label}\n→ Tag {int(day_hit)}",
                xy=(day_hit, pct_target),
                xytext=(day_hit + 0.5, pct_target - 8),
                fontsize=8,
                color="gray",
            )
    ax.set_title("Kumulative Angriffsrate (geschlossenes System)")
    ax.set_xlabel("Tag [Tage]")
    ax.set_ylabel("Kumulativ infiziert [%]")
    ax.set_ylim(0, 110)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Ensemble (mehrere Seeds)
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
    """Erstellt 4 aggregierte Visualisierungen über alle Ensemble-Läufe."""
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
        f"Einzelgitter-Kalibrierung (Ensemble)  |  {n_runs} Läufe  |"
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
    ax.set_xlabel("Transmissionsrate β₀ [pro Kontakt*Carrier]")
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
    ax.set_ylabel("Läufe [Anzahl]")
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
    ax.set_xlabel("Neue Fälle an Tag 1 [Patienten]")
    ax.set_ylabel("Anteil Läufe [%]")
    ax.set_xticks(range(int(day1_cases.max()) + 2))
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, axis="y")

    plt.tight_layout()
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _run_ensemble(
    raw: Dict[str, Any],
    config_path: Path,
    n_runs: int,
    output_dir: Path | None,
    beta_0: float,
    contact_attempts: float,
    hygiene: float,
    beta_eff: float,
    lambda_0: float,
    n_total: int,
    n_sus_0: int,
    n_car_0: int,
    days: int,
) -> None:
    """Führt n_runs Läufe mit Seeds 0..n_runs-1 aus und aggregiert die Ergebnisse."""
    out = output_dir or (
        PROJECT_ROOT
        / "outputs"
        / (datetime.now().strftime("%Y%m%d_%H%M%S") + "_Single_Ward_Ensemble")
    )
    data_dir = out / "data"
    plot_dir = out / "plots"

    lambda_target_lo, lambda_target_hi = 0.0046, 0.0054
    status = (
        "OK  -- innerhalb Zielbereich"
        if lambda_target_lo <= lambda_0 <= lambda_target_hi
        else "WARNUNG -- ausserhalb Zielbereich"
    )
    print("=" * 65)
    print(f"Einzelgitter-Kalibrierung (Ensemble): {n_runs} Läufe")
    print(f"  Config:  {config_path}")
    print(f"  β₀={beta_0:.3f}  c={contact_attempts:.1f} Kontakte/Tag  H={hygiene:.2f}")
    print(
        f"  → β_eff = {beta_eff:.4f} Tag⁻¹  |  λ(0) = {lambda_0:.4f} Tag⁻¹"
        f"  |  Ziel: {lambda_target_lo}–{lambda_target_hi} Tag⁻¹  |  {status}"
    )
    print(f"  Personen: N={n_total} (S₀={n_sus_0}, I₀={n_car_0})  |  Tage/Lauf: {days}")
    print(f"  Seeds: 0 bis {n_runs - 1} (deterministisch, reproduzierbar)")
    print("=" * 65)

    all_dfs: List[pd.DataFrame] = []
    all_metas: List[Dict[str, Any]] = []
    for i in range(n_runs):
        df, meta = run_once(raw, seed=i)
        all_dfs.append(df.assign(run_id=i))
        all_metas.append(meta)
        if (i + 1) % max(1, n_runs // 4) == 0 or i == 0:
            pct_done = (i + 1) / n_runs * 100
            print(f"  [{pct_done:5.1f}%] Lauf {i + 1:>{len(str(n_runs))}}/{n_runs} abgeschlossen")

    print("=" * 65)
    print("Aggregiere Ergebnisse...")
    combined = pd.concat(all_dfs, ignore_index=True)
    summaries = pd.DataFrame(all_metas)

    for milestone in [25, 50, 75, 100]:
        vals = summaries[f"days_to_{milestone}pct"].dropna()
        if len(vals) > 0:
            print(
                f"  {milestone:3d}% erreicht in {len(vals)}/{n_runs} Läufen"
                f"  |  Median: {vals.median():.0f} Tage"
                f"  |  [P5={vals.quantile(0.05):.0f}, P95={vals.quantile(0.95):.0f}]"
            )
        else:
            print(f"  {milestone:3d}% nicht erreicht in {days} Tagen (0/{n_runs} Läufe)")

    data_dir.mkdir(parents=True, exist_ok=True)
    combined_path = data_dir / "ensemble_all_runs.parquet"
    summary_path = data_dir / "ensemble_summary.parquet"
    combined.to_parquet(combined_path, index=False)
    summaries.to_parquet(summary_path, index=False)
    print("\nGespeichert:")
    print(f"  Alle Läufe → {combined_path}")
    print(f"  Zusammenfassung → {summary_path}")

    plot_path = plot_dir / "01_single_ward_ensemble.png"
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


# ---------------------------------------------------------------------------
# Einstiegspunkt (Einzellauf oder Ensemble via --n-runs)
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Einzelgitter-Kalibrierung für Transmissionsparameter"
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
        default=1,
        help="Anzahl Läufe. 1 = detaillierter Einzellauf, >1 = Ensemble über Seeds 0..N-1",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Ausgabeverzeichnis (default: outputs/YYYYMMDD_HHMMSS)",
    )
    args = parser.parse_args()

    raw = _load_config(args.config)

    run_cfg = raw.get("run", {})
    days: int = int(run_cfg.get("days", 100))
    seed: int = int(run_cfg.get("seed", 42))

    pop_cfg = raw.get("population", {})
    n_sus_0: int = int(pop_cfg.get("susceptible_count", 8))
    n_car_0: int = int(pop_cfg.get("carrier_count", 2))
    n_total: int = n_sus_0 + n_car_0

    beta_0, contact_attempts, hygiene, beta_eff, _ = _derive_transmission_params(raw)
    lambda_0 = beta_eff * n_car_0 / n_total

    if args.n_runs > 1:
        _run_ensemble(
            raw=raw,
            config_path=args.config,
            n_runs=args.n_runs,
            output_dir=args.output_dir,
            beta_0=beta_0,
            contact_attempts=contact_attempts,
            hygiene=hygiene,
            beta_eff=beta_eff,
            lambda_0=lambda_0,
            n_total=n_total,
            n_sus_0=n_sus_0,
            n_car_0=n_car_0,
            days=days,
        )
        return

    output_dir = args.output_dir or (
        PROJECT_ROOT / "outputs" / (datetime.now().strftime("%Y%m%d_%H%M%S") + "_Single_Ward")
    )
    data_dir = output_dir / "data"
    plot_dir = output_dir / "plots"

    print("=" * 60)
    print("Einzelgitter-Kalibrierung")
    print(f"  Config:  {args.config}")
    print(f"  β₀={beta_0:.3f}  c={contact_attempts:.1f} Kontakte/Tag  H={hygiene:.2f}")
    lambda_target_lo, lambda_target_hi = 0.0046, 0.0054
    in_range = lambda_target_lo <= lambda_0 <= lambda_target_hi
    status = "OK  -- innerhalb Zielbereich" if in_range else "WARNUNG -- ausserhalb Zielbereich"
    print(f"  → β_eff = {beta_eff:.4f} Tag⁻¹")
    print(
        f"  → λ(0)  = {lambda_0:.4f} Tag⁻¹  |  Ziel: {lambda_target_lo}–{lambda_target_hi} Tag⁻¹  |  {status}"
    )
    print(f"  Personen: N={n_total} (S₀={n_sus_0}, I₀={n_car_0})")
    print(f"  Seed: {seed}  |  Tage: {days}")
    print("=" * 60)

    df, meta = run_once(raw, seed)

    # Tabelle ausgeben (bis 100% oder Ende)
    stop_day = meta["days_to_100pct"] or days
    print(f"{'Tag':>5} {'S(t)':>6} {'I(t)':>6} {'Neu':>5} {'I/N [%]':>9} {'Kumul. [%]':>11}")
    print("-" * 50)
    for _, row in df[df["day"] <= stop_day].iterrows():
        print(
            f"{int(row.day):>5}  {int(row.susceptible):>5}  {int(row.carriers):>5}"
            f"  {int(row.new_cases):>4}  {row.prevalence_pct:>8.1f}%  {row.cumulative_attacked_pct:>10.1f}%"
        )

    print("=" * 60)
    print("\nMeilensteine (geschlossenes System):")
    for milestone in [25, 50, 75, 100]:
        day_hit = meta[f"days_to_{milestone}pct"]
        if day_hit is not None:
            print(f"  {milestone:3d}% → Tag {day_hit}")
        else:
            print(f"  {milestone:3d}% → nicht erreicht in {days} Tagen")

    data_dir.mkdir(parents=True, exist_ok=True)
    parquet_path = data_dir / "calibration_daily.parquet"
    df.to_parquet(parquet_path, index=False)
    print("\nGespeichert:")
    print(f"  Daten → {parquet_path}")

    plot_path = plot_dir / "01_single_ward_calibration.png"
    _plot(
        df=df,
        beta_eff=beta_eff,
        beta_0=beta_0,
        contact_attempts=contact_attempts,
        hygiene=hygiene,
        n_total=n_total,
        n_sus_0=n_sus_0,
        n_car_0=n_car_0,
        seed=seed,
        plot_path=plot_path,
    )
    print(f"  Plot  → {plot_path}")


if __name__ == "__main__":
    main()
