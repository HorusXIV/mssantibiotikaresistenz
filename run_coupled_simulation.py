"""Run the coupled macro + micro simulation from ``shared/config.yml``."""

from __future__ import annotations

import argparse
import csv
import math
import random
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import yaml

import visualize_results
from exchange.patient import Department, HealthState, Patient
from macro_simulation.simulation import SimulationConfig as MacroConfig
from macro_simulation.simulator import MacroSimulator
from micro_simulation.simulation import SimulationConfig as MicroConfig
from micro_simulation.simulator import MicroSimulator

DEFAULT_CONFIG_PATH = Path(__file__).resolve().parent / "shared" / "config_realistic.yml"
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / "outputs"
DEFAULT_CSV_DIR = DEFAULT_OUTPUT_DIR / "csv"
_FORBIDDEN_TEMPLATE_FIELDS = {
    "_ctx",
    "admission_day",
    "department",
    "episode_id",
    "hospital_id",
    "is_isolated",
    "patient_id",
    "planned_discharge_day",
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


_DAILY_FIELDS = [
    "day",
    "run_id",
    "total_patients",
    "susceptible",
    "carriers",
    "prevalence",
    "avg_resistant_fraction",
    "isolated_count",
    "abx_on_count",
    "ward_count",
    "icu_count",
    "isolation_count",
]

_DAILY_BY_HOSPITAL_FIELDS = [
    "day",
    "run_id",
    "hospital_id",
    "total_patients",
    "susceptible",
    "carriers",
    "prevalence",
    "avg_resistant_fraction",
    "isolated_count",
    "abx_on_count",
    "ward_count",
    "icu_count",
    "isolation_count",
]

_MICRO_DAILY_FIELDS = [
    "day",
    "run_id",
    "carrier_count",
    "active_episodes",
    "abx_on_carrier_count",
    "isolated_carrier_count",
    "mean_resistant_fraction",
    "p10_resistant_fraction",
    "p50_resistant_fraction",
    "p90_resistant_fraction",
    "mean_p_clearance",
    "mean_relative_transmissibility",
    "mean_n_strains",
    "mean_total_population",
    "genotype_entropy",
    "genotype_S_fraction",
    "genotype_R1_fraction",
    "genotype_R2_fraction",
    "genotype_R3_fraction",
    "genotype_other_fraction",
]

_MICRO_DAILY_BY_HOSPITAL_FIELDS = [
    "day",
    "run_id",
    "hospital_id",
    "carrier_count",
    "active_episodes",
    "abx_on_carrier_count",
    "isolated_carrier_count",
    "mean_resistant_fraction",
    "mean_p_clearance",
    "mean_relative_transmissibility",
    "mean_n_strains",
    "mean_total_population",
    "genotype_entropy",
]

_MICRO_PATIENT_DAILY_FIELDS = [
    "day",
    "run_id",
    "hospital_id",
    "patient_id",
    "episode_id",
    "department",
    "is_isolated",
    "abx_on",
    "abx_class",
    "dose_level",
    "adherence",
    "immune_strength",
    "immune_status",
    "resistant_fraction",
    "dominant_genotype",
    "relative_transmissibility",
    "p_clearance",
    "severity_modifier",
    "lethality_modifier",
    "episode_day",
    "n_strains",
    "total_population",
]

_MICRO_DAILY_GENOTYPE_FIELDS = [
    "day",
    "run_id",
    "dominant_genotype",
    "count",
    "fraction",
]


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


def _make_patient_factory(
    population: PopulationSettings,
    macro_config: MacroConfig,
    seed: int,
) -> Callable[[str, Department], Patient]:
    """Return a closure that creates new patients for the admission phase."""
    counter = [0]
    rng = random.Random(seed + 99999)

    def factory(_hospital_id: str, _department: Department) -> Patient:
        counter[0] += 1
        pid = f"dyn_{counter[0]:07d}"
        is_carrier = rng.random() < macro_config.community_carrier_fraction
        if is_carrier:
            return Patient(
                patient_id=pid,
                state=HealthState.CARRIER,
                episode_id=f"community_ep_{counter[0]:07d}",
                resistant_fraction=macro_config.replacement_resistant_fraction,
                dominant_genotype=macro_config.replacement_dominant_genotype,
                **population.susceptible_template,
            )
        return Patient(
            patient_id=pid,
            state=HealthState.SUSCEPTIBLE,
            **population.susceptible_template,
        )

    return factory


def _collect_patient_stats(patients: list[Patient]) -> dict[str, float | int]:
    total = len(patients)
    carriers = [p for p in patients if p.state == HealthState.CARRIER]
    susceptible = total - len(carriers)
    avg_res = sum(p.resistant_fraction for p in carriers) / len(carriers) if carriers else 0.0
    isolated_count = sum(1 for p in patients if p.is_isolated)
    abx_on_count = sum(1 for p in patients if p.regimen.on)
    ward_count = sum(1 for p in patients if p.department == Department.WARD)
    icu_count = sum(1 for p in patients if p.department == Department.ICU)
    isolation_count = sum(1 for p in patients if p.department == Department.ISOLATION)
    prevalence = len(carriers) / total if total else 0.0

    return {
        "total_patients": total,
        "susceptible": susceptible,
        "carriers": len(carriers),
        "prevalence": prevalence,
        "avg_resistant_fraction": avg_res,
        "isolated_count": isolated_count,
        "abx_on_count": abx_on_count,
        "ward_count": ward_count,
        "icu_count": icu_count,
        "isolation_count": isolation_count,
    }


def _collect_macro_daily_logs(
    macro: MacroSimulator,
    n_hospitals: int,
    day: int,
    run_id: str,
) -> tuple[dict[str, float | int | str], list[dict[str, float | int | str]]]:
    all_patients: list[Patient] = []
    by_hospital: list[dict[str, float | int | str]] = []

    for i in range(1, n_hospitals + 1):
        hid = f"hospital_{i:03d}"
        patients = macro.get_patients(hid)
        all_patients.extend(patients)
        row = {
            "day": day,
            "run_id": run_id,
            "hospital_id": hid,
        }
        row.update(_collect_patient_stats(patients))
        by_hospital.append(row)

    global_row: dict[str, float | int | str] = {
        "day": day,
        "run_id": run_id,
    }
    global_row.update(_collect_patient_stats(all_patients))

    return global_row, by_hospital


def _mean_or_zero(values: list[float]) -> float:
    if not values:
        return 0.0
    return float(sum(values) / len(values))


def _quantile_or_zero(values: list[float], q: float) -> float:
    if not values:
        return 0.0
    xs = sorted(values)
    if len(xs) == 1:
        return float(xs[0])
    pos = q * (len(xs) - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return float(xs[lo])
    weight = pos - lo
    return float(xs[lo] * (1.0 - weight) + xs[hi] * weight)


def _aggregate_micro_records(records: list[dict[str, Any]]) -> dict[str, float | int]:
    carrier_count = len(records)
    if carrier_count == 0:
        return {
            "carrier_count": 0,
            "active_episodes": 0,
            "abx_on_carrier_count": 0,
            "isolated_carrier_count": 0,
            "mean_resistant_fraction": 0.0,
            "p10_resistant_fraction": 0.0,
            "p50_resistant_fraction": 0.0,
            "p90_resistant_fraction": 0.0,
            "mean_p_clearance": 0.0,
            "mean_relative_transmissibility": 0.0,
            "mean_n_strains": 0.0,
            "mean_total_population": 0.0,
            "genotype_entropy": 0.0,
            "genotype_S_fraction": 0.0,
            "genotype_R1_fraction": 0.0,
            "genotype_R2_fraction": 0.0,
            "genotype_R3_fraction": 0.0,
            "genotype_other_fraction": 0.0,
        }

    resistant_values = [float(r["resistant_fraction"]) for r in records]
    clearance_values = [float(r["p_clearance"]) for r in records]
    transmissibility_values = [float(r["relative_transmissibility"]) for r in records]
    n_strain_values = [float(r["n_strains"]) for r in records if bool(r["has_state"])]
    total_population_values = [
        float(r["total_population"]) for r in records if bool(r["has_state"])
    ]

    genotype_counts = Counter(str(r["dominant_genotype"]) for r in records)
    fractions = {g: c / carrier_count for g, c in genotype_counts.items()}
    entropy = -sum(frac * math.log2(frac) for frac in fractions.values() if frac > 0.0)
    tracked = {"S", "R1", "R2", "R3"}
    other_fraction = sum(frac for genotype, frac in fractions.items() if genotype not in tracked)

    return {
        "carrier_count": carrier_count,
        "active_episodes": sum(1 for r in records if bool(r["has_state"])),
        "abx_on_carrier_count": sum(1 for r in records if bool(r["abx_on"])),
        "isolated_carrier_count": sum(1 for r in records if bool(r["is_isolated"])),
        "mean_resistant_fraction": _mean_or_zero(resistant_values),
        "p10_resistant_fraction": _quantile_or_zero(resistant_values, 0.10),
        "p50_resistant_fraction": _quantile_or_zero(resistant_values, 0.50),
        "p90_resistant_fraction": _quantile_or_zero(resistant_values, 0.90),
        "mean_p_clearance": _mean_or_zero(clearance_values),
        "mean_relative_transmissibility": _mean_or_zero(transmissibility_values),
        "mean_n_strains": _mean_or_zero(n_strain_values),
        "mean_total_population": _mean_or_zero(total_population_values),
        "genotype_entropy": float(entropy),
        "genotype_S_fraction": float(fractions.get("S", 0.0)),
        "genotype_R1_fraction": float(fractions.get("R1", 0.0)),
        "genotype_R2_fraction": float(fractions.get("R2", 0.0)),
        "genotype_R3_fraction": float(fractions.get("R3", 0.0)),
        "genotype_other_fraction": float(other_fraction),
    }


def _collect_micro_daily_logs(
    macro: MacroSimulator,
    micro: MicroSimulator,
    n_hospitals: int,
    day: int,
    run_id: str,
) -> tuple[
    dict[str, float | int | str],
    list[dict[str, float | int | str]],
    list[dict[str, float | int | str]],
    list[dict[str, float | int | str]],
]:
    all_records: list[dict[str, Any]] = []
    by_hospital_rows: list[dict[str, float | int | str]] = []
    patient_rows: list[dict[str, float | int | str]] = []
    genotype_counter: Counter[str] = Counter()

    for i in range(1, n_hospitals + 1):
        hid = f"hospital_{i:03d}"
        carriers = [p for p in macro.get_patients(hid) if p.state == HealthState.CARRIER]
        hospital_records: list[dict[str, Any]] = []

        for patient in carriers:
            episode_day = 0
            n_strains = 0
            total_population = 0.0
            has_state = False

            episode_id = patient.episode_id
            if episode_id:
                episode_state = micro.get_episode_state(episode_id)
                if episode_state is not None:
                    has_state = True
                    episode_day = int(episode_state.day)
                    n_strains = int(episode_state.population.n_strains)
                    total_population = float(episode_state.population.total_population)

            record = {
                "dominant_genotype": patient.dominant_genotype,
                "resistant_fraction": float(patient.resistant_fraction),
                "p_clearance": float(patient.p_clearance),
                "relative_transmissibility": float(patient.relative_transmissibility),
                "n_strains": n_strains,
                "total_population": total_population,
                "has_state": has_state,
                "abx_on": bool(patient.regimen.on),
                "is_isolated": bool(patient.is_isolated),
            }
            hospital_records.append(record)
            all_records.append(record)
            genotype_counter[str(patient.dominant_genotype)] += 1

            patient_rows.append(
                {
                    "day": day,
                    "run_id": run_id,
                    "hospital_id": hid,
                    "patient_id": patient.patient_id,
                    "episode_id": episode_id or "",
                    "department": patient.department.value,
                    "is_isolated": bool(patient.is_isolated),
                    "abx_on": bool(patient.regimen.on),
                    "abx_class": patient.regimen.abx_class,
                    "dose_level": patient.regimen.dose_level,
                    "adherence": float(patient.adherence),
                    "immune_strength": float(patient.immune_strength),
                    "immune_status": patient.immune_status,
                    "resistant_fraction": float(patient.resistant_fraction),
                    "dominant_genotype": patient.dominant_genotype,
                    "relative_transmissibility": float(patient.relative_transmissibility),
                    "p_clearance": float(patient.p_clearance),
                    "severity_modifier": float(patient.severity_modifier),
                    "lethality_modifier": float(patient.lethality_modifier),
                    "episode_day": episode_day,
                    "n_strains": n_strains,
                    "total_population": total_population,
                }
            )

        hospital_row = {
            "day": day,
            "run_id": run_id,
            "hospital_id": hid,
        }
        stats = _aggregate_micro_records(hospital_records)
        hospital_row.update(
            {
                "carrier_count": stats["carrier_count"],
                "active_episodes": stats["active_episodes"],
                "abx_on_carrier_count": stats["abx_on_carrier_count"],
                "isolated_carrier_count": stats["isolated_carrier_count"],
                "mean_resistant_fraction": stats["mean_resistant_fraction"],
                "mean_p_clearance": stats["mean_p_clearance"],
                "mean_relative_transmissibility": stats["mean_relative_transmissibility"],
                "mean_n_strains": stats["mean_n_strains"],
                "mean_total_population": stats["mean_total_population"],
                "genotype_entropy": stats["genotype_entropy"],
            }
        )
        by_hospital_rows.append(hospital_row)

    global_row: dict[str, float | int | str] = {"day": day, "run_id": run_id}
    global_stats = _aggregate_micro_records(all_records)
    global_row.update(global_stats)
    global_row["active_episodes"] = len(micro.get_active_episodes())

    total_carriers = int(global_stats["carrier_count"])
    genotype_rows: list[dict[str, float | int | str]] = []
    genotype_order = ["S", "R1", "R2", "R3"]
    for genotype in genotype_order:
        count = int(genotype_counter.get(genotype, 0))
        genotype_rows.append(
            {
                "day": day,
                "run_id": run_id,
                "dominant_genotype": genotype,
                "count": count,
                "fraction": (count / total_carriers) if total_carriers > 0 else 0.0,
            }
        )
    other_count = sum(
        count for genotype, count in genotype_counter.items() if genotype not in genotype_order
    )
    genotype_rows.append(
        {
            "day": day,
            "run_id": run_id,
            "dominant_genotype": "OTHER",
            "count": other_count,
            "fraction": (other_count / total_carriers) if total_carriers > 0 else 0.0,
        }
    )

    return global_row, by_hospital_rows, patient_rows, genotype_rows


def _write_csv(path: Path, rows: list[dict[str, float | int | str]], fieldnames: list[str]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


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

    patient_factory = _make_patient_factory(
        population=settings.population,
        macro_config=settings.macro,
        seed=settings.run.seed,
    )

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
    macro_daily_rows: list[dict[str, float | int | str]] = []
    macro_daily_by_hospital_rows: list[dict[str, float | int | str]] = []
    micro_daily_rows: list[dict[str, float | int | str]] = []
    micro_daily_by_hospital_rows: list[dict[str, float | int | str]] = []
    micro_patient_daily_rows: list[dict[str, float | int | str]] = []
    micro_daily_genotype_rows: list[dict[str, float | int | str]] = []
    for day in range(1, settings.run.days + 1):
        macro.step(
            micro_simulator=micro,
            run_id=settings.run.run_id,
            patient_factory=patient_factory,
        )
        summary = _summarize_day(
            macro=macro,
            n_hospitals=settings.population.hospitals,
            day=day,
        )
        final_summary = summary

        global_row, hospital_rows = _collect_macro_daily_logs(
            macro=macro,
            n_hospitals=settings.population.hospitals,
            day=day,
            run_id=settings.run.run_id,
        )
        macro_daily_rows.append(global_row)
        macro_daily_by_hospital_rows.extend(hospital_rows)
        (
            micro_global_row,
            micro_hospital_rows,
            micro_pat_rows,
            micro_genotype_rows,
        ) = _collect_micro_daily_logs(
            macro=macro,
            micro=micro,
            n_hospitals=settings.population.hospitals,
            day=day,
            run_id=settings.run.run_id,
        )
        micro_daily_rows.append(micro_global_row)
        micro_daily_by_hospital_rows.extend(micro_hospital_rows)
        micro_patient_daily_rows.extend(micro_pat_rows)
        micro_daily_genotype_rows.extend(micro_genotype_rows)

        if not settings.run.quiet:
            print(
                f"day={summary.day:03d} susceptible={summary.susceptible:04d} "
                f"carriers={summary.carriers:04d} "
                f"avg_resistant_fraction={summary.avg_resistant_fraction:.4f}"
            )

    if final_summary is None:
        return

    macro_daily_path = DEFAULT_CSV_DIR / "macro_daily.csv"
    macro_daily_by_hospital_path = DEFAULT_CSV_DIR / "macro_daily_by_hospital.csv"
    micro_daily_path = DEFAULT_CSV_DIR / "micro_daily.csv"
    micro_daily_by_hospital_path = DEFAULT_CSV_DIR / "micro_daily_by_hospital.csv"
    micro_patient_daily_path = DEFAULT_CSV_DIR / "micro_patient_daily.csv"
    micro_daily_genotype_path = DEFAULT_CSV_DIR / "micro_daily_genotype.csv"
    _write_csv(macro_daily_path, macro_daily_rows, _DAILY_FIELDS)
    _write_csv(
        macro_daily_by_hospital_path,
        macro_daily_by_hospital_rows,
        _DAILY_BY_HOSPITAL_FIELDS,
    )
    _write_csv(micro_daily_path, micro_daily_rows, _MICRO_DAILY_FIELDS)
    _write_csv(
        micro_daily_by_hospital_path,
        micro_daily_by_hospital_rows,
        _MICRO_DAILY_BY_HOSPITAL_FIELDS,
    )
    _write_csv(micro_patient_daily_path, micro_patient_daily_rows, _MICRO_PATIENT_DAILY_FIELDS)
    _write_csv(micro_daily_genotype_path, micro_daily_genotype_rows, _MICRO_DAILY_GENOTYPE_FIELDS)

    print(
        "run_end "
        f"day={final_summary.day} susceptible={final_summary.susceptible} "
        f"carriers={final_summary.carriers} "
        f"avg_resistant_fraction={final_summary.avg_resistant_fraction:.4f} "
        f"active_micro_episodes={len(micro.get_active_episodes())}"
    )
    print(f"macro_log_written daily={macro_daily_path} by_hospital={macro_daily_by_hospital_path}")
    print(
        "micro_log_written "
        f"daily={micro_daily_path} "
        f"by_hospital={micro_daily_by_hospital_path} "
        f"patient_daily={micro_patient_daily_path} "
        f"genotype={micro_daily_genotype_path}"
    )

    visualize_results.run(
        csv_dir=DEFAULT_CSV_DIR,
        plot_dir=DEFAULT_OUTPUT_DIR / "plots",
        quiet=True,
    )


if __name__ == "__main__":
    main()
