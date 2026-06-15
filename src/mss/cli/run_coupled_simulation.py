"""Run the coupled macro + micro simulation from ``config/simulation_realistic.yml``."""

from __future__ import annotations

import argparse
import copy
import importlib.metadata
import json
import math
import platform
import random
import subprocess
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import yaml

from mss.cli import visualize_results
from mss.domain import Department, HealthState, Patient
from mss.simulation.macro import MacroSimulator
from mss.simulation.macro import SimulationConfig as MacroConfig
from mss.simulation.micro import GeneIndex, MicroSimulator, build_micro_config, classify_genotype
from mss.simulation.micro import SimulationConfig as MicroConfig

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CONFIG_PATH = PROJECT_ROOT / "config" / "simulation_realistic.yml"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "outputs"
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
    use_micro: bool = True  # when False, mss-run skips the micro layer entirely


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
]

_MACRO_CELL_DAILY_FIELDS = [
    "day",
    "run_id",
    "hospital_id",
    "x",
    "y",
    "department",
    "total_patients",
    "carriers",
    "prevalence",
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
    "resistant_fraction",
    "dominant_genotype",
    "dominant_strain_name",
    "relative_transmissibility",
    "p_clearance",
    "severity_modifier",
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

_MICRO_HOSPITAL_POPULATION_FIELDS = [
    "day",
    "snapshot_stage",
    "run_id",
    "hospital_id",
    "carrier_count",
    "tracked_episode_count",
    "total_population",
    "total_strains",
]

_MICRO_STRAIN_DAILY_FIELDS = [
    "day",
    "snapshot_stage",
    "run_id",
    "hospital_id",
    "patient_id",
    "episode_id",
    "episode_day",
    "strain_id",
    "parent_id",
    "donor_id",
    "founder_id",
    "strain_name",
    "strain_genotype",
    "is_dominant",
    "population",
    "population_fraction_in_episode",
    "lineage_age",
    "damage_load",
]

_MICRO_EPISODE_GENE_DAILY_FIELDS = [
    "day",
    "snapshot_stage",
    "run_id",
    "hospital_id",
    "patient_id",
    "episode_id",
    "episode_day",
    "gene_index",
    "gene_name",
    "weighted_mean",
    "weighted_variance",
    "present_strain_count",
    "total_strain_count",
    "present_strain_fraction",
    "present_population_fraction",
    "dominant_strain_value",
    "presence_threshold",
]
_TRANSFER_DAILY_FIELDS = [
    "day",
    "run_id",
    "from_hospital",
    "to_hospital",
    "transfers",
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
        seed=_require_non_negative_int(run_raw.get("seed"), "run.seed"),
        run_id=str(run_raw.get("run_id", "dev_run")),
        quiet=bool(run_raw.get("quiet", False)),
        use_micro=bool(run_raw.get("use_micro", True)),
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
    micro = build_micro_config(micro_raw, source="micro")

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
            # Community carriers use the carrier template (calibrated immunity,
            # compliance, severity, clearance); only the resistance state is
            # overridden with the community-import values.
            template = dict(population.carrier_template)
            template["resistant_fraction"] = macro_config.replacement_resistant_fraction
            template["dominant_genotype"] = macro_config.replacement_dominant_genotype
            return Patient(
                patient_id=pid,
                state=HealthState.CARRIER,
                episode_id=f"community_ep_{counter[0]:07d}",
                **template,
            )
        return Patient(
            patient_id=pid,
            state=HealthState.SUSCEPTIBLE,
            **population.susceptible_template,
        )

    return factory


def _iter_all_patients(macro: MacroSimulator, n_hospitals: int) -> list[tuple[str, Patient]]:
    rows: list[tuple[str, Patient]] = []
    for i in range(1, n_hospitals + 1):
        hid = f"hospital_{i:03d}"
        rows.extend((hid, patient) for patient in macro.get_patients(hid))
    return rows


def _seed_initial_micro_states(
    macro: MacroSimulator,
    micro: MicroSimulator,
    n_hospitals: int,
    seed: int,
) -> None:
    """Materialize day-0 micro states for already-admitted carrier episodes."""
    rng = random.Random(seed + 424242)
    carriers = sorted(
        (
            (hid, patient)
            for hid, patient in _iter_all_patients(macro, n_hospitals)
            if patient.state == HealthState.CARRIER and patient.episode_id
        ),
        key=lambda item: (item[0], item[1].patient_id),
    )
    for _hid, patient in carriers:
        micro.initialize_episode(
            episode_id=str(patient.episode_id),
            patient_id=patient.patient_id,
            resistant_fraction=float(patient.resistant_fraction),
            dominant_genotype=str(patient.dominant_genotype),
            dominant_strain_name=patient.dominant_strain_name or None,
            day=0,
            seed=rng.randint(0, 2**31 - 1),
        )


def _collect_patient_stats(patients: list[Patient]) -> dict[str, float | int]:
    total = len(patients)
    carriers = [p for p in patients if p.state == HealthState.CARRIER]
    susceptible = total - len(carriers)
    avg_res = sum(p.resistant_fraction for p in carriers) / len(carriers) if carriers else 0.0
    isolated_count = sum(1 for p in patients if p.is_isolated)
    abx_on_count = sum(1 for p in patients if p.regimen.on)
    ward_count = sum(1 for p in patients if p.department == Department.WARD)
    icu_count = sum(1 for p in patients if p.department == Department.ICU)
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


def _collect_cell_daily_logs(
    macro: MacroSimulator,
    n_hospitals: int,
    day: int,
    run_id: str,
) -> list[dict[str, float | int | str]]:
    """Per-cell occupancy rows for every hospital grid cell (carrier prevalence)."""
    rows: list[dict[str, float | int | str]] = []
    for i in range(1, n_hospitals + 1):
        hid = f"hospital_{i:03d}"
        for cell in macro.get_cell_stats(hid):
            rows.append(
                {
                    "day": day,
                    "run_id": run_id,
                    "hospital_id": hid,
                    "x": int(cell["x"]),
                    "y": int(cell["y"]),
                    "department": str(cell["department"]),
                    "total_patients": int(cell["total_patients"]),
                    "carriers": int(cell["carriers"]),
                    "prevalence": float(cell["prevalence"]),
                }
            )
    return rows


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
                    "resistant_fraction": float(patient.resistant_fraction),
                    "dominant_genotype": patient.dominant_genotype,
                    "dominant_strain_name": patient.dominant_strain_name,
                    "relative_transmissibility": float(patient.relative_transmissibility),
                    "p_clearance": float(patient.p_clearance),
                    "severity_modifier": float(patient.severity_modifier),
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


def _collect_micro_raw_logs(
    macro: MacroSimulator,
    micro: MicroSimulator,
    n_hospitals: int,
    day: int,
    run_id: str,
    snapshot_stage: str,
) -> tuple[
    list[dict[str, float | int | str]],
    list[dict[str, float | int | str]],
    list[dict[str, float | int | str]],
]:
    hospital_population_rows: list[dict[str, float | int | str]] = []
    strain_rows: list[dict[str, float | int | str]] = []
    gene_rows: list[dict[str, float | int | str]] = []
    presence_threshold = float(micro.config.gene_presence_threshold)

    for i in range(1, n_hospitals + 1):
        hid = f"hospital_{i:03d}"
        carriers = [
            patient for patient in macro.get_patients(hid) if patient.state == HealthState.CARRIER
        ]
        tracked_episode_count = 0
        total_population = 0.0
        total_strains = 0

        for patient in carriers:
            if not patient.episode_id:
                continue
            episode_state = micro.get_episode_state(patient.episode_id)
            if episode_state is None:
                continue

            tracked_episode_count += 1
            population = episode_state.population
            total_population += float(population.total_population)
            total_strains += int(population.n_strains)

            if population.n_strains == 0 or population.total_population <= 0.0:
                continue

            dominant_idx = int(np.argmax(population.populations))
            pop_total = float(population.total_population)
            for idx in range(population.n_strains):
                strain_population = float(population.populations[idx])
                if strain_population <= 0.0:
                    continue
                strain_rows.append(
                    {
                        "day": day,
                        "snapshot_stage": snapshot_stage,
                        "run_id": run_id,
                        "hospital_id": hid,
                        "patient_id": patient.patient_id,
                        "episode_id": patient.episode_id,
                        "episode_day": int(episode_state.day),
                        "strain_id": population.strain_ids[idx],
                        "parent_id": population.parent_ids[idx],
                        "donor_id": population.donor_ids[idx],
                        "founder_id": population.founder_ids[idx],
                        "strain_name": population.strain_names[idx],
                        "strain_genotype": classify_genotype(population.genomes[idx]),
                        "is_dominant": idx == dominant_idx,
                        "population": strain_population,
                        "population_fraction_in_episode": strain_population / pop_total,
                        "lineage_age": float(population.lineage_ages[idx]),
                        "damage_load": float(population.damage_loads[idx]),
                    }
                )

            weights = population.populations.astype(np.float64, copy=False)
            genomes = population.genomes.astype(np.float64, copy=False)
            dominant_genome = population.genomes[dominant_idx]
            for gene_idx in range(len(GeneIndex)):
                gene_values = genomes[:, gene_idx]
                weighted_mean = float(np.average(gene_values, weights=weights))
                weighted_variance = float(
                    np.average(np.square(gene_values - weighted_mean), weights=weights)
                )
                present_mask = gene_values >= presence_threshold
                present_strain_count = int(np.sum(present_mask))
                present_population_fraction = float(np.sum(weights[present_mask]) / pop_total)
                gene_rows.append(
                    {
                        "day": day,
                        "snapshot_stage": snapshot_stage,
                        "run_id": run_id,
                        "hospital_id": hid,
                        "patient_id": patient.patient_id,
                        "episode_id": patient.episode_id,
                        "episode_day": int(episode_state.day),
                        "gene_index": gene_idx,
                        "gene_name": GeneIndex(gene_idx).name.lower(),
                        "weighted_mean": weighted_mean,
                        "weighted_variance": weighted_variance,
                        "present_strain_count": present_strain_count,
                        "total_strain_count": int(population.n_strains),
                        "present_strain_fraction": (
                            present_strain_count / float(population.n_strains)
                            if population.n_strains > 0
                            else 0.0
                        ),
                        "present_population_fraction": present_population_fraction,
                        "dominant_strain_value": float(dominant_genome[gene_idx]),
                        "presence_threshold": presence_threshold,
                    }
                )

        hospital_population_rows.append(
            {
                "day": day,
                "snapshot_stage": snapshot_stage,
                "run_id": run_id,
                "hospital_id": hid,
                "carrier_count": len(carriers),
                "tracked_episode_count": tracked_episode_count,
                "total_population": total_population,
                "total_strains": total_strains,
            }
        )

    return hospital_population_rows, strain_rows, gene_rows


def _write_parquet(
    path: Path, rows: list[dict[str, float | int | str]], fieldnames: list[str]
) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows, columns=fieldnames)
    df.to_parquet(path, index=False)


class _ParquetBatchWriter:
    """Append parquet row groups incrementally to keep memory bounded."""

    def __init__(self, path: Path, fieldnames: list[str]) -> None:
        self.path = path
        self.fieldnames = fieldnames
        self._writer: pq.ParquetWriter | None = None

    def append_rows(self, rows: list[dict[str, float | int | str]]) -> None:
        if not rows:
            return

        self.path.parent.mkdir(parents=True, exist_ok=True)
        df = pd.DataFrame(rows, columns=self.fieldnames)
        table = pa.Table.from_pandas(df, preserve_index=False)
        if self._writer is None:
            self._writer = pq.ParquetWriter(self.path, table.schema)
        self._writer.write_table(table)

    def close(self) -> None:
        if self._writer is not None:
            self._writer.close()
            self._writer = None


def _build_settings_from_raw(raw: dict[str, Any], seed: int) -> CoupledSimulationSettings:
    """Build CoupledSimulationSettings from an in-memory config dict.

    Mirrors load_coupled_settings() but accepts a dict instead of a file path.
    Always sets quiet=True and overrides the seed — intended for calibration loops.
    """
    raw = copy.deepcopy(raw)
    run_raw = _require_mapping("run", raw.get("run"))
    population_raw = _require_mapping("population", raw.get("population"))
    macro_raw = dict(_require_mapping("macro", raw.get("macro")))
    micro_raw = dict(_require_mapping("micro", raw.get("micro")))

    micro_workers = micro_raw.pop("workers", None)
    if micro_workers is not None:
        micro_workers = _require_positive_int(micro_workers, "micro.workers")

    run = RunSettings(
        days=_require_positive_int(run_raw.get("days"), "run.days"),
        seed=seed,
        run_id=str(run_raw.get("run_id", "sweep")),
        quiet=True,
    )
    population = PopulationSettings(
        hospitals=_require_positive_int(population_raw.get("hospitals"), "population.hospitals"),
        susceptible_count=_require_non_negative_int(
            population_raw.get("susceptible_count"), "population.susceptible_count"
        ),
        carrier_count=_require_non_negative_int(
            population_raw.get("carrier_count"), "population.carrier_count"
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
    return CoupledSimulationSettings(
        config_path=Path("(in-memory)"),
        run=run,
        population=population,
        macro=MacroConfig(**macro_raw),
        micro=build_micro_config(micro_raw, source="micro"),
        micro_workers=micro_workers,
    )


def run_realistic_once(
    raw: dict[str, Any], seed: int, use_micro: bool = True
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Run one coupled simulation from a config dict (no files, plots, or prints).

    Analogous to run_once() in run_single_ward_calibration.py; intended for
    calibration and sweep loops.

    Args:
        use_micro: If False, the micro simulator is skipped and patients keep their
            template defaults (p_clearance, relative_transmissibility,
            severity_modifier). Recommended for the macro calibrations.

    Returns:
        df: daily DataFrame with day, total_patients, susceptible, carriers,
            prevalence, isolated_count, new_cases.
        meta: summary over the last third (seed, mean_prevalence,
            acquisition_rate_per_1000, isolation_fraction).
    """
    settings = _build_settings_from_raw(raw, seed)

    macro = MacroSimulator(
        config=settings.macro,
        n_hospitals=settings.population.hospitals,
        seed=settings.run.seed,
    )
    micro: MicroSimulator | None = (
        MicroSimulator(config=settings.micro, n_workers=1) if use_micro else None
    )

    _admit_initial_population(macro=macro, population=settings.population)
    if micro is not None:
        _seed_initial_micro_states(
            macro=macro,
            micro=micro,
            n_hospitals=settings.population.hospitals,
            seed=settings.run.seed,
        )
    patient_factory = _make_patient_factory(
        population=settings.population,
        macro_config=settings.macro,
        seed=settings.run.seed,
    )

    rows = []
    # True in-hospital acquisitions: per-day delta of the colonization counter
    # (only S->C transmission, excludes community-imported carriers and churn).
    prev_colonizations = macro.get_colonization_count()
    for day in range(1, settings.run.days + 1):
        macro.step(
            micro_simulator=micro,
            run_id=settings.run.run_id,
            patient_factory=patient_factory,
        )
        all_patients = [p for _, p in _iter_all_patients(macro, settings.population.hospitals)]
        stats = _collect_patient_stats(all_patients)
        colonizations = macro.get_colonization_count()
        new_cases = colonizations - prev_colonizations
        prev_colonizations = colonizations
        rows.append({"day": day, "new_cases": new_cases, **stats})

    if micro is not None:
        micro.close()

    df = pd.DataFrame(rows)
    last_third_start = max(1, int(settings.run.days * 2 / 3))
    last = df[df["day"] >= last_third_start]

    mean_carriers = float(last["carriers"].mean())
    mean_susceptible = float(last["susceptible"].mean())
    mean_new_cases = float(last["new_cases"].mean())
    mean_isolated = float(last["isolated_count"].mean())

    acquisition_rate = mean_new_cases / mean_susceptible * 1000 if mean_susceptible > 0 else 0.0
    isolation_fraction = mean_isolated / mean_carriers if mean_carriers > 0 else 0.0

    # Spatial transmission structure: fraction of in-hospital transmissions whose
    # source was a same-cell (roommate) carrier — set by proximity_decay_alpha.
    n_transmissions = macro.get_transmission_count()
    same_cell_transmission_fraction = (
        macro.get_same_cell_transmission_count() / n_transmissions if n_transmissions > 0 else 0.0
    )

    meta: dict[str, Any] = {
        "seed": seed,
        "mean_prevalence": float(last["prevalence"].mean()),
        "acquisition_rate_per_1000": acquisition_rate,
        "isolation_fraction": isolation_fraction,
        "mean_resistant_fraction": float(last["avg_resistant_fraction"].mean()),
        "same_cell_transmission_fraction": same_cell_transmission_fraction,
        "n_transmissions": n_transmissions,
    }
    return df, meta


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
            "Defaults to config/simulation_realistic.yml in the project root."
        ),
    )
    return parser.parse_args()


def _git_provenance(root: Path) -> dict[str, Any]:
    """Best-effort git commit + dirty flag; ``None`` entries if unavailable."""

    def _run(args: list[str]) -> str | None:
        try:
            result = subprocess.run(args, cwd=root, capture_output=True, text=True, check=True)
            return result.stdout.strip()
        except (subprocess.SubprocessError, OSError):
            return None

    commit = _run(["git", "rev-parse", "HEAD"])
    status = _run(["git", "status", "--porcelain"])
    return {"commit": commit, "dirty": None if status is None else bool(status)}


def _json_default(obj: Any) -> Any:
    """Fallback JSON encoder for numpy/enum/path/set values in the run metadata."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, (set, frozenset)):
        return sorted(obj)
    if isinstance(obj, Department):
        return obj.value
    if isinstance(obj, Path):
        return str(obj)
    return str(obj)


def _write_run_meta(
    data_dir: Path,
    settings: CoupledSimulationSettings,
    run_ts: str,
) -> Path:
    """Persist a provenance record so a run is traceable to its exact inputs."""
    try:
        package_version = importlib.metadata.version("mss")
    except importlib.metadata.PackageNotFoundError:
        package_version = "unknown"

    meta = {
        "run_timestamp": run_ts,
        "config_path": str(settings.config_path),
        "git": _git_provenance(PROJECT_ROOT),
        "python_version": platform.python_version(),
        "package_version": package_version,
        "seed": settings.run.seed,
        "run": asdict(settings.run),
        "population": asdict(settings.population),
        "macro": asdict(settings.macro),
        "micro": asdict(settings.micro),
        "micro_workers": settings.micro_workers,
    }
    path = data_dir / "run_meta.json"
    with path.open("w", encoding="utf-8") as handle:
        json.dump(meta, handle, indent=2, default=_json_default, ensure_ascii=False)
    return path


def run(
    settings: CoupledSimulationSettings,
    output_dir: Path | None = None,
) -> Path:
    """Run the coupled simulation with the given settings and write outputs.

    Args:
        settings: The complete simulation setup.
        output_dir: If provided, use this directory for outputs. If None, a
            timestamped directory is created in DEFAULT_OUTPUT_DIR.

    Returns:
        The path to the run directory.
    """
    run_ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = output_dir or (DEFAULT_OUTPUT_DIR / run_ts)
    data_dir = run_dir / "data"
    plot_dir = run_dir / "plots"
    data_dir.mkdir(parents=True, exist_ok=True)

    run_meta_path = _write_run_meta(data_dir, settings, run_ts)
    if not settings.run.quiet:
        print(f"run_meta {run_meta_path}")

    macro = MacroSimulator(
        config=settings.macro,
        n_hospitals=settings.population.hospitals,
        seed=settings.run.seed,
    )
    micro = MicroSimulator(
        config=settings.micro,
        n_workers=settings.micro_workers,
    )

    # When micro is disabled, no episodes are materialized and macro.step receives
    # micro_simulator=None; carriers keep their template defaults throughout the run.
    coupled_micro = micro if settings.run.use_micro else None

    _admit_initial_population(macro=macro, population=settings.population)
    if coupled_micro is not None:
        _seed_initial_micro_states(
            macro=macro,
            micro=micro,
            n_hospitals=settings.population.hospitals,
            seed=settings.run.seed,
        )

    patient_factory = _make_patient_factory(
        population=settings.population,
        macro_config=settings.macro,
        seed=settings.run.seed,
    )

    if not settings.run.quiet:
        print(
            "run_start "
            f"config={settings.config_path} days={settings.run.days} "
            f"hospitals={settings.population.hospitals} "
            f"susceptible={settings.population.susceptible_count} "
            f"seed_carriers={settings.population.carrier_count} "
            f"micro={'on' if settings.run.use_micro else 'off'} "
            f"micro_steps_per_day={settings.micro.steps_per_day} "
            f"micro_workers={micro.n_workers} seed={settings.run.seed}"
        )

    final_summary = None
    macro_daily_path = data_dir / "macro_daily.parquet"
    macro_daily_by_hospital_path = data_dir / "macro_daily_by_hospital.parquet"
    macro_cell_daily_path = data_dir / "macro_cell_daily.parquet"
    micro_daily_path = data_dir / "micro_daily.parquet"
    micro_daily_by_hospital_path = data_dir / "micro_daily_by_hospital.parquet"
    micro_patient_daily_path = data_dir / "micro_patient_daily.parquet"
    micro_daily_genotype_path = data_dir / "micro_daily_genotype.parquet"
    micro_hospital_population_path = data_dir / "micro_hospital_population.parquet"
    micro_strain_daily_path = data_dir / "micro_strain_daily.parquet"
    micro_episode_gene_daily_path = data_dir / "micro_episode_gene_daily.parquet"
    transfer_daily_path = data_dir / "transfer_daily.parquet"

    macro_daily_writer = _ParquetBatchWriter(macro_daily_path, _DAILY_FIELDS)
    macro_daily_by_hospital_writer = _ParquetBatchWriter(
        macro_daily_by_hospital_path,
        _DAILY_BY_HOSPITAL_FIELDS,
    )
    macro_cell_daily_writer = _ParquetBatchWriter(
        macro_cell_daily_path,
        _MACRO_CELL_DAILY_FIELDS,
    )
    micro_daily_writer = _ParquetBatchWriter(micro_daily_path, _MICRO_DAILY_FIELDS)
    micro_daily_by_hospital_writer = _ParquetBatchWriter(
        micro_daily_by_hospital_path,
        _MICRO_DAILY_BY_HOSPITAL_FIELDS,
    )
    micro_patient_daily_writer = _ParquetBatchWriter(
        micro_patient_daily_path,
        _MICRO_PATIENT_DAILY_FIELDS,
    )
    micro_daily_genotype_writer = _ParquetBatchWriter(
        micro_daily_genotype_path,
        _MICRO_DAILY_GENOTYPE_FIELDS,
    )
    micro_hospital_population_writer = _ParquetBatchWriter(
        micro_hospital_population_path,
        _MICRO_HOSPITAL_POPULATION_FIELDS,
    )
    micro_strain_daily_writer = _ParquetBatchWriter(
        micro_strain_daily_path,
        _MICRO_STRAIN_DAILY_FIELDS,
    )
    micro_episode_gene_daily_writer = _ParquetBatchWriter(
        micro_episode_gene_daily_path,
        _MICRO_EPISODE_GENE_DAILY_FIELDS,
    )
    transfer_daily_writer = _ParquetBatchWriter(transfer_daily_path, _TRANSFER_DAILY_FIELDS)

    (
        initial_hospital_population_rows,
        initial_strain_rows,
        initial_gene_rows,
    ) = _collect_micro_raw_logs(
        macro=macro,
        micro=micro,
        n_hospitals=settings.population.hospitals,
        day=0,
        run_id=settings.run.run_id,
        snapshot_stage="initial",
    )
    micro_hospital_population_writer.append_rows(initial_hospital_population_rows)
    micro_strain_daily_writer.append_rows(initial_strain_rows)
    micro_episode_gene_daily_writer.append_rows(initial_gene_rows)

    for day in range(1, settings.run.days + 1):
        macro.step(
            micro_simulator=coupled_micro,
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
        macro_daily_writer.append_rows([global_row])
        macro_daily_by_hospital_writer.append_rows(hospital_rows)
        macro_cell_daily_writer.append_rows(
            _collect_cell_daily_logs(
                macro=macro,
                n_hospitals=settings.population.hospitals,
                day=day,
                run_id=settings.run.run_id,
            )
        )
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
        micro_daily_writer.append_rows([micro_global_row])
        micro_daily_by_hospital_writer.append_rows(micro_hospital_rows)
        micro_patient_daily_writer.append_rows(micro_pat_rows)
        micro_daily_genotype_writer.append_rows(micro_genotype_rows)
        (
            hospital_population_rows,
            strain_rows,
            gene_rows,
        ) = _collect_micro_raw_logs(
            macro=macro,
            micro=micro,
            n_hospitals=settings.population.hospitals,
            day=day,
            run_id=settings.run.run_id,
            snapshot_stage="end_of_day",
        )
        micro_hospital_population_writer.append_rows(hospital_population_rows)
        micro_strain_daily_writer.append_rows(strain_rows)
        micro_episode_gene_daily_writer.append_rows(gene_rows)

        transfer_counts = Counter(
            (t["from_hospital"], t["to_hospital"]) for t in macro.get_daily_transfers()
        )
        transfer_rows = []
        for (src, dst), n in transfer_counts.items():
            transfer_rows.append(
                {
                    "day": day,
                    "run_id": settings.run.run_id,
                    "from_hospital": src,
                    "to_hospital": dst,
                    "transfers": n,
                }
            )
        transfer_daily_writer.append_rows(transfer_rows)

        if not settings.run.quiet:
            print(
                f"day={summary.day:03d} susceptible={summary.susceptible:04d} "
                f"carriers={summary.carriers:04d} "
                f"avg_resistant_fraction={summary.avg_resistant_fraction:.4f}"
            )

    macro_daily_writer.close()
    macro_daily_by_hospital_writer.close()
    macro_cell_daily_writer.close()
    micro_daily_writer.close()
    micro_daily_by_hospital_writer.close()
    micro_patient_daily_writer.close()
    micro_daily_genotype_writer.close()
    micro_hospital_population_writer.close()
    micro_strain_daily_writer.close()
    micro_episode_gene_daily_writer.close()
    transfer_daily_writer.close()

    if final_summary is not None:
        if not settings.run.quiet:
            print(
                "run_end "
                f"day={final_summary.day} susceptible={final_summary.susceptible} "
                f"carriers={final_summary.carriers} "
                f"avg_resistant_fraction={final_summary.avg_resistant_fraction:.4f} "
                f"active_micro_episodes={len(micro.get_active_episodes())}"
            )
            print(f"run_dir={run_dir}")
            print(
                f"macro_log_written daily={macro_daily_path} by_hospital={macro_daily_by_hospital_path}"
            )
            print(
                "micro_log_written "
                f"daily={micro_daily_path} "
                f"by_hospital={micro_daily_by_hospital_path} "
                f"patient_daily={micro_patient_daily_path} "
                f"genotype={micro_daily_genotype_path} "
                f"hospital_population={micro_hospital_population_path} "
                f"strain_daily={micro_strain_daily_path} "
                f"episode_gene_daily={micro_episode_gene_daily_path}"
            )

    visualize_results.run(
        data_dir=data_dir,
        plot_dir=plot_dir,
        quiet=True,
    )

    return run_dir


def main() -> None:
    """Load the config, run the coupled simulation, and write outputs and plots."""
    args = _parse_args()
    config_path = args.config if args.config is not None else DEFAULT_CONFIG_PATH
    settings = load_coupled_settings(config_path)
    run(settings)


if __name__ == "__main__":
    main()
