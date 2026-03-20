"""CLI runner for coupled macro + micro simulation.

Usage example:
    uv run python run_coupled_simulation.py --days 60 --hospitals 3 --seed 42
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import List

from exchange.patient import Department, HealthState, Patient
from macro_simulation.simulation import SimulationConfig as MacroConfig
from macro_simulation.simulator import MacroSimulator
from micro_simulation.simulation import SimulationConfig as MicroConfig
from micro_simulation.simulator import MicroSimulator


@dataclass
class DaySummary:
    day: int
    susceptible: int
    carriers: int
    avg_resistant_fraction: float


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run coupled macro+micro MSS simulation.")

    parser.add_argument("--days", type=int, default=60, help="Number of macro days to simulate.")
    parser.add_argument("--hospitals", type=int, default=1, help="Number of hospitals in network.")
    parser.add_argument(
        "--susceptible",
        type=int,
        default=100,
        help="Number of initially susceptible patients.",
    )
    parser.add_argument(
        "--seed-carriers",
        type=int,
        default=5,
        help="Number of initially carrier patients.",
    )
    parser.add_argument(
        "--steps-per-day",
        type=int,
        default=12,
        help="Micro steps per macro day (default 12).",
    )
    parser.add_argument("--seed", type=int, default=42, help="RNG seed for reproducible runs.")
    parser.add_argument(
        "--run-id",
        type=str,
        default="dev_run",
        help="Run identifier passed to micro requests.",
    )
    parser.add_argument(
        "--department",
        choices=["ward", "icu"],
        default="ward",
        help="Initial department for all admitted patients.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress per-day logs and print only final summary.",
    )

    return parser


def _hospital_id_for(i: int, n_hospitals: int) -> str:
    return f"hospital_{(i % n_hospitals) + 1:03d}"


def _admit_initial_population(
    macro: MacroSimulator,
    n_hospitals: int,
    n_susceptible: int,
    n_seed_carriers: int,
    department: Department,
) -> None:
    for i in range(n_susceptible):
        patient = Patient(patient_id=f"sus_{i:05d}", state=HealthState.SUSCEPTIBLE)
        macro.admit(
            patient,
            hospital_id=_hospital_id_for(i, n_hospitals),
            department=department,
        )

    for i in range(n_seed_carriers):
        patient = Patient(
            patient_id=f"car_{i:05d}",
            state=HealthState.CARRIER,
            episode_id=f"seed_ep_{i:05d}",
            resistant_fraction=0.2,
            dominant_genotype="R1",
            relative_transmissibility=1.2,
            p_clearance=0.03,
        )
        macro.admit(
            patient,
            hospital_id=_hospital_id_for(i + n_susceptible, n_hospitals),
            department=department,
        )


def _summarize_day(macro: MacroSimulator, n_hospitals: int, day: int) -> DaySummary:
    patients: List[Patient] = []
    for i in range(1, n_hospitals + 1):
        patients.extend(macro.get_patients(f"hospital_{i:03d}"))

    carriers = [p for p in patients if p.state == HealthState.CARRIER]
    susceptible = len(patients) - len(carriers)
    avg_res = sum(p.resistant_fraction for p in carriers) / len(carriers) if carriers else 0.0

    return DaySummary(
        day=day,
        susceptible=susceptible,
        carriers=len(carriers),
        avg_resistant_fraction=avg_res,
    )


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    if args.days <= 0:
        raise ValueError("--days must be > 0")
    if args.hospitals <= 0:
        raise ValueError("--hospitals must be > 0")
    if args.susceptible < 0:
        raise ValueError("--susceptible must be >= 0")
    if args.seed_carriers < 0:
        raise ValueError("--seed-carriers must be >= 0")
    if args.steps_per_day <= 0:
        raise ValueError("--steps-per-day must be > 0")

    department = Department.WARD if args.department == "ward" else Department.ICU

    macro = MacroSimulator(config=MacroConfig(), n_hospitals=args.hospitals, seed=args.seed)
    micro = MicroSimulator(config=MicroConfig(steps_per_day=args.steps_per_day))

    _admit_initial_population(
        macro=macro,
        n_hospitals=args.hospitals,
        n_susceptible=args.susceptible,
        n_seed_carriers=args.seed_carriers,
        department=department,
    )

    print(
        "run_start "
        f"days={args.days} hospitals={args.hospitals} susceptible={args.susceptible} "
        f"seed_carriers={args.seed_carriers} micro_steps_per_day={args.steps_per_day} seed={args.seed}"
    )

    final_summary = None
    for day in range(1, args.days + 1):
        macro.step(micro_simulator=micro, run_id=args.run_id)
        summary = _summarize_day(macro=macro, n_hospitals=args.hospitals, day=day)
        final_summary = summary

        if not args.quiet:
            print(
                f"day={summary.day:03d} susceptible={summary.susceptible:04d} "
                f"carriers={summary.carriers:04d} avg_resistant_fraction={summary.avg_resistant_fraction:.4f}"
            )

    if final_summary is None:
        return

    print(
        "run_end "
        f"day={final_summary.day} susceptible={final_summary.susceptible} "
        f"carriers={final_summary.carriers} "
        f"avg_resistant_fraction={final_summary.avg_resistant_fraction:.4f} "
        f"active_micro_episodes={len(micro.get_active_episodes())}"
    )


if __name__ == "__main__":
    main()
