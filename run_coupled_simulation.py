"""Run the coupled macro + micro simulation from ``shared/config.yml``."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from exchange.patient import Department, HealthState, Patient
from macro_simulation.simulation import SimulationConfig as MacroConfig
from macro_simulation.simulator import MacroSimulator
from micro_simulation.simulation import SimulationConfig as MicroConfig
from micro_simulation.simulator import MicroSimulator

DEFAULT_CONFIG_PATH = Path(__file__).resolve().parent / "shared" / "config_realistic.yml"
_FORBIDDEN_TEMPLATE_FIELDS = {
    "_ctx",
    "department",
    "episode_id",
    "hospital_id",
    "is_isolated",
    "patient_id",
    "regimen",
    "state",
}


@dataclass
class DaySummary:
    day: int
    susceptible: int
    carriers: int
    avg_resistant_fraction: float


@dataclass
class RunSettings:
    days: int
    seed: int
    run_id: str
    quiet: bool


@dataclass
class PopulationSettings:
    hospitals: int
    susceptible_count: int
    carrier_count: int
    initial_department: Department
    susceptible_template: dict[str, Any]
    carrier_template: dict[str, Any]


@dataclass
class CoupledSimulationSettings:
    config_path: Path
    run: RunSettings
    population: PopulationSettings
    macro: MacroConfig
    micro: MicroConfig
    micro_workers: int | None


def _require_mapping(section_name: str, value: Any) -> dict[str, Any]:
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise TypeError(f"Config section '{section_name}' must be a mapping.")
    return dict(value)


def _require_positive_int(value: Any, field_name: str) -> int:
    if not isinstance(value, int) or value <= 0:
        raise ValueError(f"Config field '{field_name}' must be a positive integer.")
    return value


def _require_non_negative_int(value: Any, field_name: str) -> int:
    if not isinstance(value, int) or value < 0:
        raise ValueError(f"Config field '{field_name}' must be a non-negative integer.")
    return value


def _parse_department(value: Any) -> Department:
    if value == Department.WARD.value:
        return Department.WARD
    if value == Department.ICU.value:
        return Department.ICU
    raise ValueError("Config field 'population.initial_department' must be 'ward' or 'icu'.")


def _normalize_patient_template(section_name: str, template: dict[str, Any]) -> dict[str, Any]:
    forbidden = sorted(_FORBIDDEN_TEMPLATE_FIELDS.intersection(template))
    if forbidden:
        joined = ", ".join(forbidden)
        raise ValueError(
            f"Config section '{section_name}' contains runner-owned patient fields: {joined}."
        )

    normalized = dict(template)
    if "history_flags" in normalized:
        flags = normalized["history_flags"]
        if not isinstance(flags, list):
            raise TypeError(f"Config field '{section_name}.history_flags' must be a list.")
        normalized["history_flags"] = set(flags)

    return normalized


def load_coupled_settings(config_path: Path = DEFAULT_CONFIG_PATH) -> CoupledSimulationSettings:
    """Load the entire coupled simulation setup from YAML."""
    if not config_path.exists():
        raise FileNotFoundError(f"Simulation config not found: {config_path}")

    raw = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict):
        raise TypeError("Top-level config must be a mapping.")

    run_raw = _require_mapping("run", raw.get("run"))
    population_raw = _require_mapping("population", raw.get("population"))
    macro_raw = _require_mapping("macro", raw.get("macro"))
    micro_raw = _require_mapping("micro", raw.get("micro"))

    run = RunSettings(
        days=_require_positive_int(run_raw.get("days"), "run.days"),
        seed=_require_positive_int(run_raw.get("seed"), "run.seed"),
        run_id=str(run_raw.get("run_id", "dev_run")),
        quiet=bool(run_raw.get("quiet", False)),
    )

    micro_workers = micro_raw.pop("workers", None)
    if micro_workers is not None:
        micro_workers = _require_positive_int(micro_workers, "micro.workers")

    population = PopulationSettings(
        hospitals=_require_positive_int(population_raw.get("hospitals"), "population.hospitals"),
        susceptible_count=_require_non_negative_int(
            population_raw.get("susceptible_count"),
            "population.susceptible_count",
        ),
        carrier_count=_require_non_negative_int(
            population_raw.get("carrier_count"),
            "population.carrier_count",
        ),
        initial_department=_parse_department(population_raw.get("initial_department")),
        susceptible_template=_normalize_patient_template(
            "population.susceptible_template",
            _require_mapping(
                "population.susceptible_template",
                population_raw.get("susceptible_template"),
            ),
        ),
        carrier_template=_normalize_patient_template(
            "population.carrier_template",
            _require_mapping("population.carrier_template", population_raw.get("carrier_template")),
        ),
    )

    macro = MacroConfig(**macro_raw)
    micro = MicroConfig(**micro_raw)

    return CoupledSimulationSettings(
        config_path=config_path,
        run=run,
        population=population,
        macro=macro,
        micro=micro,
        micro_workers=micro_workers,
    )


def _hospital_id_for(i: int, n_hospitals: int) -> str:
    return f"hospital_{(i % n_hospitals) + 1:03d}"


def _build_patient(
    patient_id: str,
    state: HealthState,
    template: dict[str, Any],
    episode_id: str | None = None,
) -> Patient:
    patient = Patient(
        patient_id=patient_id,
        state=state,
        episode_id=episode_id,
        **template,
    )
    return patient


def _admit_initial_population(
    macro: MacroSimulator,
    population: PopulationSettings,
) -> None:
    for i in range(population.susceptible_count):
        patient = _build_patient(
            patient_id=f"sus_{i:05d}",
            state=HealthState.SUSCEPTIBLE,
            template=population.susceptible_template,
        )
        macro.admit(
            patient,
            hospital_id=_hospital_id_for(i, population.hospitals),
            department=population.initial_department,
        )

    for i in range(population.carrier_count):
        patient = _build_patient(
            patient_id=f"car_{i:05d}",
            state=HealthState.CARRIER,
            episode_id=f"seed_ep_{i:05d}",
            template=population.carrier_template,
        )
        macro.admit(
            patient,
            hospital_id=_hospital_id_for(i + population.susceptible_count, population.hospitals),
            department=population.initial_department,
        )


def _summarize_day(macro: MacroSimulator, n_hospitals: int, day: int) -> DaySummary:
    patients: list[Patient] = []
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


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the coupled macro + micro antibiotic resistance simulation.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        metavar="FILE",
        help=(
            "Path to the YAML config file. "
            "Defaults to shared/config_realistic.yml next to this script."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    config_path = args.config if args.config is not None else DEFAULT_CONFIG_PATH

    settings = load_coupled_settings(config_path)

    macro = MacroSimulator(
        config=settings.macro,
        n_hospitals=settings.population.hospitals,
        seed=settings.run.seed,
    )
    micro = MicroSimulator(
        config=settings.micro,
        n_workers=settings.micro_workers,
    )

    _admit_initial_population(macro=macro, population=settings.population)

    print(
        "run_start "
        f"config={settings.config_path} days={settings.run.days} "
        f"hospitals={settings.population.hospitals} "
        f"susceptible={settings.population.susceptible_count} "
        f"seed_carriers={settings.population.carrier_count} "
        f"micro_steps_per_day={settings.micro.steps_per_day} "
        f"micro_workers={micro.n_workers} seed={settings.run.seed}"
    )

    final_summary = None
    for day in range(1, settings.run.days + 1):
        macro.step(micro_simulator=micro, run_id=settings.run.run_id)
        summary = _summarize_day(
            macro=macro,
            n_hospitals=settings.population.hospitals,
            day=day,
        )
        final_summary = summary

        if not settings.run.quiet:
            print(
                f"day={summary.day:03d} susceptible={summary.susceptible:04d} "
                f"carriers={summary.carriers:04d} "
                f"avg_resistant_fraction={summary.avg_resistant_fraction:.4f}"
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
