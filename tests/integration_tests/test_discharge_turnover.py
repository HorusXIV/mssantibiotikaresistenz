from __future__ import annotations

import math
import random

import pytest

from mss.domain import Department, HealthState, Patient
from mss.simulation.macro import MacroSimulator, SimulationConfig

_H1 = "hospital_001"


class _RecordingMicro:
    def __init__(self) -> None:
        self.cleared: list[str] = []

    def clear_episode(self, episode_id: str) -> None:
        self.cleared.append(episode_id)


class _FixedRng:
    def __init__(self, value: float) -> None:
        self._value = value

    def random(self) -> float:
        return self._value


def test_discharge_fix_micro_clears_episode() -> None:
    sim = MacroSimulator(config=SimulationConfig(), n_hospitals=1, seed=1)
    patient = Patient(
        patient_id="car_001",
        state=HealthState.CARRIER,
        episode_id="ep_001",
    )
    sim.admit(patient, hospital_id=_H1, department=Department.WARD)

    micro = _RecordingMicro()
    sim.discharge(patient, micro_simulator=micro)

    assert micro.cleared == ["ep_001"]


def test_discharge_preserves_carrier_state() -> None:
    sim = MacroSimulator(config=SimulationConfig(), n_hospitals=1, seed=2)
    patient = Patient(
        patient_id="car_002",
        state=HealthState.CARRIER,
        episode_id="ep_002",
        resistant_fraction=0.4,
        dominant_genotype="R1",
    )
    sim.admit(patient, hospital_id=_H1, department=Department.WARD)
    sim.discharge(patient)

    assert patient.state == HealthState.CARRIER
    assert patient.episode_id == "ep_002"
    assert patient.resistant_fraction == pytest.approx(0.4)
    assert patient.dominant_genotype == "R1"


def test_draw_los_ward_mean_approx_correct() -> None:
    cfg = SimulationConfig(los_mean_ward=8.0, los_sigma=0.6)
    sim = MacroSimulator(config=cfg, n_hospitals=1, seed=3)

    samples = [sim._draw_los(Department.WARD) for _ in range(4000)]
    mean = sum(samples) / len(samples)

    assert 7.0 <= mean <= 9.5


def test_draw_los_icu_longer_than_ward() -> None:
    cfg = SimulationConfig(los_mean_ward=8.0, los_mean_icu=14.0, los_sigma=0.6)
    sim = MacroSimulator(config=cfg, n_hospitals=1, seed=4)

    ward_mean = sum(sim._draw_los(Department.WARD) for _ in range(2000)) / 2000
    icu_mean = sum(sim._draw_los(Department.ICU) for _ in range(2000)) / 2000

    assert icu_mean > ward_mean


def test_carrier_extends_planned_discharge_day() -> None:
    cfg = SimulationConfig(daily_admission_rate=1.0, carrier_extension_days=10)
    sim = MacroSimulator(config=cfg, n_hospitals=1, seed=5)
    patient = Patient(
        patient_id="car_ext",
        state=HealthState.CARRIER,
        episode_id="ep_ext",
    )
    sim.admit(patient, hospital_id=_H1, department=Department.WARD)
    sim._day = 7
    patient.planned_discharge_day = 7

    sim._discharge_phase([patient], _H1, micro_simulator=None)

    assert patient.planned_discharge_day == 17


def test_logistic_p_increases_over_target() -> None:
    cfg = SimulationConfig(discharge_logistic_k=1.0, discharge_logistic_t_half=3.0)

    def logistic(days_over: int) -> float:
        return 1.0 / (
            1.0 + math.exp(-cfg.discharge_logistic_k * (days_over - cfg.discharge_logistic_t_half))
        )

    assert logistic(0) == pytest.approx(0.0474, rel=0.2)
    assert logistic(3) == pytest.approx(0.5, rel=0.05)
    assert logistic(7) >= 0.9


def test_susceptible_after_target_eventually_discharged() -> None:
    cfg = SimulationConfig(daily_admission_rate=1.0)
    sim = MacroSimulator(config=cfg, n_hospitals=1, seed=6)
    patient = Patient(patient_id="sus_ext", state=HealthState.SUSCEPTIBLE)
    sim.admit(patient, hospital_id=_H1, department=Department.WARD)

    sim._rng = _FixedRng(0.6)
    patient.planned_discharge_day = 0

    discharged = False
    for day in range(10):
        sim._day = day
        patients = sim.get_patients(_H1)
        if not patients:
            discharged = True
            break
        sim._discharge_phase(patients, _H1, micro_simulator=None)
        discharged = patient.hospital_id is None
        if discharged:
            break

    assert discharged


def test_admission_phase_poisson_mean_respected() -> None:
    rng = random.Random(7)
    lam = 8.0
    draws = [MacroSimulator._poisson_draw(rng, lam) for _ in range(5000)]
    mean = sum(draws) / len(draws)

    assert 7.2 <= mean <= 8.8
