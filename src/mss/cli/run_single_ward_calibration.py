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
    mss-calibrate --config config/simulation_single_ward.yml

Aufruf (1000 Laeufe):
    mss-calibrate-batch --n-runs 1000
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CONFIG = PROJECT_ROOT / "config" / "simulation_single_ward.yml"


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

    # --- [0,0] SI-Dynamik ---
    ax = axes[0, 0]
    ax.plot(days, df["susceptible"], color="steelblue", lw=2, label="Susceptible S(t)")
    ax.plot(days, carriers, color="firebrick", lw=2, label="Carrier I(t)")
    ax.axhline(n_total, color="gray", lw=1, ls="--", alpha=0.5, label=f"Total N={n_total}")
    ax.set_title("Transmission Dynamics")
    ax.set_xlabel("Tag [Tage]")
    ax.set_ylabel("Patienten")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    ax.set_ylim(0, n_total * 1.15)

    # --- [0,1] Taeglich Inzidenz ---
    ax = axes[0, 1]
    ax.bar(days, new_cases, color="firebrick", alpha=0.7, label="Neue Faelle/Tag")
    ax.set_title("Tägliche Inzidenz")
    ax.set_xlabel("Tag [Tage]")
    ax.set_ylabel("Neue Fälle [Tag⁻¹]")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, axis="y")
    if n_total > 0:
        ax2 = ax.twinx()
        ax2.set_ylabel(f"pro 100 Patienten-Tage (N={n_total})", fontsize=8, color="gray")
        ax2.set_ylim(0, max(ax.get_ylim()[1] / n_total * 100, 0.01))
        ax2.tick_params(colors="gray")

    # --- [1,0] Infektionskraft lambda(t) ---
    ax = axes[1, 0]
    ax.plot(days, df["lambda_t"], color="darkorange", lw=2, label="λ(t) simuliert")
    ax.axhline(
        lambda_0_theoretical,
        color="darkorange",
        lw=1,
        ls="--",
        alpha=0.6,
        label=f"λ(0) theoretisch = {lambda_0_theoretical:.4f} Tag⁻¹",
    )
    ax.set_title("Force of Infection λ(t)")
    ax.set_xlabel("Tag [Tage]")
    ax.set_ylabel("Infektionshazard [Tag⁻¹]")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # --- [1,1] Kumulative Angriffsrate ---
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
    ax.set_title("Cumulative Attack Rate")
    ax.set_xlabel("Tag [Tage]")
    ax.set_ylabel("Kumulativ infiziert [%]")
    ax.set_ylim(0, 110)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot gespeichert → {plot_path}")


# ---------------------------------------------------------------------------
# Einzellauf-Einstiegspunkt
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Einzelgitter-Kalibrierung fuer Transmissionsparameter (Einzellauf)"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="Pfad zur YAML-Konfigurationsdatei",
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

    output_dir = args.output_dir or (
        PROJECT_ROOT / "outputs" / datetime.now().strftime("%Y%m%d_%H%M%S")
    )
    data_dir = output_dir / "data"
    plot_dir = output_dir / "plots"

    print("=" * 60)
    print("Einzelgitter-Kalibrierung")
    print(f"  Config:  {args.config}")
    print(f"  β₀={beta_0:.3f}  c={contact_attempts:.1f} Kontakte/Tag  H={hygiene:.2f}")
    print(f"  → β_eff = {beta_eff:.4f} Tag⁻¹")
    print(f"  → λ(0)  = {lambda_0:.4f} Tag⁻¹")
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
    if meta["days_to_100pct"] is not None:
        print(f"→ 100% Infektion erreicht an Tag {meta['days_to_100pct']}")
    else:
        print(
            f"→ Simulation beendet nach {days} Tagen ({meta['final_attack_pct']:.1f}% kumulativ infiziert)"
        )

    print("\nMeilensteine:")
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
