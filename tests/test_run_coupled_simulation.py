from __future__ import annotations

import json
from pathlib import Path

import yaml

from mss.cli.run_coupled_simulation import (
    _admit_initial_population,
    _collect_micro_raw_logs,
    _seed_initial_micro_states,
    _write_run_meta,
    load_coupled_settings,
)
from mss.domain import Department, HealthState
from mss.simulation.macro import MacroSimulator
from mss.simulation.micro import MicroSimulator


def _write_config(path: Path) -> None:
    config = {
        "run": {
            "days": 3,
            "seed": 7,
            "run_id": "test_run",
            "quiet": True,
        },
        "population": {
            "hospitals": 2,
            "susceptible_count": 4,
            "carrier_count": 2,
            "initial_department": "icu",
            "susceptible_template": {
                "age_years": 33,
                "history_flags": ["prior_abx"],
            },
            "carrier_template": {
                "age_years": 70,
                "immune_strength": 0.5,
                "resistant_fraction": 0.6,
                "dominant_genotype": "R2",
                "relative_transmissibility": 1.9,
                "p_clearance": 0.01,
                "history_flags": ["recent_surgery"],
            },
        },
        "macro": {
            "base_hygiene": 0.8,
            "base_isolation_effectiveness": 0.7,
            "base_diagnostic_speed": 0.4,
            "base_transmission_rate": 0.03,
            "daily_contact_attempts": 60.0,
            "icu_abx_probability": 0.5,
            "ward_abx_probability": 0.1,
            "carrier_isolation_probability": 0.25,
        },
        "micro": {
            "workers": 1,
            "founder_pool_size": 12,
            "founder_pool_seed": 3,
            "founder_pool_gene_noise_std": 0.01,
            "gene_presence_threshold": 0.25,
            "steps_per_day": 10,
            "max_strains": 25,
            "carrying_capacity": 1e8,
            "min_population": 500.0,
            "clearance_threshold": 50.0,
            "base_mutation_rate": 0.02,
            "mutation_std": 0.04,
            "stress_mutation_boost": 2.0,
            "base_hgt_rate": 0.01,
            "hgt_gene_transfer_prob": 0.2,
            "selection_strength": 1.5,
            "growth_rate_per_step": 0.2,
            "death_rate_per_step": 0.08,
            "strain_prune_threshold": 2.0,
            "base_damage_per_step": 0.02,
            "replication_damage_factor": 0.15,
            "stress_damage_factor": 0.1,
            "repair_rate_per_step": 0.06,
            "age_mortality_scale": 0.01,
            "damage_mortality_scale": 0.12,
            "lifecycle_half_life_steps": 24.0,
            "max_damage_load": 4.0,
            "dormancy_growth_penalty": 0.2,
            "synergy_repair_dormancy_bonus": 0.15,
            "synergy_stress_tolerance_bonus": 0.1,
            "stochastic_threshold": 1e4,
            "stochastic_noise_scale": 0.8,
        },
    }
    path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")


def test_load_coupled_settings_reads_yaml(tmp_path: Path):
    config_path = tmp_path / "config.yml"
    _write_config(config_path)

    settings = load_coupled_settings(config_path)

    assert settings.run.days == 3
    assert settings.run.seed == 7
    assert settings.population.hospitals == 2
    assert settings.population.initial_department == Department.ICU
    assert settings.population.susceptible_template["age_years"] == 33
    assert settings.population.susceptible_template["history_flags"] == {"prior_abx"}
    assert settings.population.carrier_template["dominant_genotype"] == "R2"
    assert settings.micro.steps_per_day == 10
    assert settings.micro.founder_pool_size == 12
    assert settings.micro.gene_presence_threshold == 0.25
    assert settings.micro_workers == 1


def test_write_run_meta_records_provenance(tmp_path: Path):
    config_path = tmp_path / "config.yml"
    _write_config(config_path)
    settings = load_coupled_settings(config_path)

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    meta_path = _write_run_meta(data_dir, settings, run_ts="20260101_000000")

    assert meta_path == data_dir / "run_meta.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8"))

    assert meta["run_timestamp"] == "20260101_000000"
    assert meta["seed"] == 7
    assert meta["run"]["days"] == 3
    assert meta["micro_workers"] == 1
    assert meta["macro"]["base_isolation_effectiveness"] == 0.7
    assert meta["micro"]["steps_per_day"] == 10
    assert meta["population"]["hospitals"] == 2
    assert set(meta["git"]) == {"commit", "dirty"}
    assert "python_version" in meta


def test_use_micro_defaults_true_and_parses_false(tmp_path: Path):
    config_path = tmp_path / "config.yml"
    _write_config(config_path)

    # Omitted in _write_config, so it must default to True.
    assert load_coupled_settings(config_path).run.use_micro is True

    raw = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    raw["run"]["use_micro"] = False
    config_path.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
    assert load_coupled_settings(config_path).run.use_micro is False


def test_use_micro_false_skips_micro_coupling(tmp_path: Path):
    config_path = tmp_path / "config.yml"
    _write_config(config_path)
    settings = load_coupled_settings(config_path)
    macro = MacroSimulator(config=settings.macro, n_hospitals=settings.population.hospitals, seed=7)
    micro = MicroSimulator(config=settings.micro, n_workers=1)

    _admit_initial_population(macro, settings.population)

    # Mirror main()'s use_micro=False path: no seeding, step receives no micro.
    for _ in range(settings.run.days):
        macro.step(micro_simulator=None, run_id=settings.run.run_id, patient_factory=None)

    # No within-host episode is ever created when micro is disabled.
    assert micro.get_active_episodes() == []

    # Seeded carriers (age 70 from the carrier template) keep their template defaults
    # because the micro layer never overwrote them.
    seeded_carriers = [
        p
        for hid in ("hospital_001", "hospital_002")
        for p in macro.get_patients(hid)
        if p.state == HealthState.CARRIER and p.age_years == 70
    ]
    assert seeded_carriers
    assert all(p.p_clearance == 0.01 for p in seeded_carriers)
    assert all(p.relative_transmissibility == 1.9 for p in seeded_carriers)


def test_admit_initial_population_uses_configured_templates(tmp_path: Path):
    config_path = tmp_path / "config.yml"
    _write_config(config_path)
    settings = load_coupled_settings(config_path)
    macro = MacroSimulator(config=settings.macro, n_hospitals=settings.population.hospitals, seed=7)

    _admit_initial_population(macro, settings.population)

    all_patients = []
    for hid in ("hospital_001", "hospital_002"):
        all_patients.extend(macro.get_patients(hid))

    carriers = [p for p in all_patients if p.state == HealthState.CARRIER]
    susceptible = [p for p in all_patients if p.state == HealthState.SUSCEPTIBLE]

    assert len(carriers) == 2
    assert len(susceptible) == 4
    assert all(p.department == Department.ICU for p in all_patients)
    assert all(p.age_years == 33 for p in susceptible)
    assert all(p.age_years == 70 for p in carriers)
    assert all(p.dominant_genotype == "R2" for p in carriers)
    assert all(p.resistant_fraction == 0.6 for p in carriers)


def test_seed_initial_micro_states_and_collect_raw_logs(tmp_path: Path):
    config_path = tmp_path / "config.yml"
    _write_config(config_path)
    settings = load_coupled_settings(config_path)
    macro = MacroSimulator(config=settings.macro, n_hospitals=settings.population.hospitals, seed=7)
    micro = MicroSimulator(config=settings.micro, n_workers=1)

    _admit_initial_population(macro, settings.population)
    _seed_initial_micro_states(
        macro=macro,
        micro=micro,
        n_hospitals=settings.population.hospitals,
        seed=settings.run.seed,
    )

    hospital_rows, strain_rows, gene_rows = _collect_micro_raw_logs(
        macro=macro,
        micro=micro,
        n_hospitals=settings.population.hospitals,
        day=0,
        run_id=settings.run.run_id,
        snapshot_stage="initial",
    )

    assert len(hospital_rows) == settings.population.hospitals
    assert (
        sum(int(row["carrier_count"]) for row in hospital_rows) == settings.population.carrier_count
    )
    assert any(float(row["total_population"]) > 0.0 for row in hospital_rows)
    assert strain_rows, "Expected raw per-strain rows for seeded carrier episodes."
    assert gene_rows, "Expected per-episode gene rows for seeded carrier episodes."
    assert all(str(row["strain_id"]).startswith("seed_ep_") for row in strain_rows)
    assert any(str(row["founder_id"]).startswith("founder_") for row in strain_rows)
    assert {row["gene_name"] for row in gene_rows}
