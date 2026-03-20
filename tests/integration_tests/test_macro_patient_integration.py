"""
Integration tests for Patient <-> MacroSimulator (hospital/network) interface compatibility.

These tests verify that:
1. The macro layer uses Patient as its sole interface object for patient state
2. Hospital placement, admissions, discharges, and transfers work correctly
3. Within-hospital transmission respects patient and hospital modifiers
4. Hygiene, isolation, antibiotic policies, and clearance are handled correctly
5. Multi-day and multi-hospital (network) simulations are stable and reproducible
6. The macro layer never depends on micro internals — only on Patient fields/methods

The macro layer is a *carrier model*: patient states are SUSCEPTIBLE and CARRIER only.
Macro owns S→C and C→S transitions. Micro-derived attributes stored on Patient
(relative_transmissibility, p_clearance, severity_modifier, lethality_modifier,
resistant_fraction) may influence macro behaviour but are never written by macro
except through Patient.clear_carriage().
"""

from __future__ import annotations

import random
from typing import Callable

import pytest

from exchange.patient import (
    AntibioticRegimen,
    Department,
    HealthState,
    Patient,
    PatientDailyContext,
)
from macro_simulation.simulation import SimulationConfig
from macro_simulation.simulator import MacroSimulator

# -- valid states in the macro carrier model ----------------------------------
_VALID_STATES = (HealthState.SUSCEPTIBLE, HealthState.CARRIER)

# -- hospital IDs used throughout tests ---------------------------------------
_HOSPITAL_IDS = [f"hospital_{i:03d}" for i in range(1, 11)]
_H1, _H2, _H3, _H4 = _HOSPITAL_IDS[:4]


# =============================================================================
# Helpers
# =============================================================================


def _make_patients(
    n_susceptible: int,
    n_carrier: int,
    *,
    carrier_p_clearance: float = 0.05,
) -> list[Patient]:
    """Create a mixed batch of susceptible and carrier patients."""
    patients: list[Patient] = []
    for i in range(n_susceptible):
        patients.append(Patient(patient_id=f"sus_{i:04d}", state=HealthState.SUSCEPTIBLE))
    for i in range(n_carrier):
        patients.append(
            Patient(
                patient_id=f"car_{i:04d}",
                state=HealthState.CARRIER,
                episode_id=f"episode_car_{i:04d}",
                resistant_fraction=0.3,
                dominant_genotype="R1",
                relative_transmissibility=1.2,
                p_clearance=carrier_p_clearance,
            )
        )
    return patients


def _count_new_carriers_mc(
    *,
    n_runs: int,
    n_days: int,
    n_susceptible: int,
    make_carrier: Callable[[int], Patient],
    config_factory: Callable[[], SimulationConfig] = SimulationConfig,
    seed_offset: int = 0,
) -> int:
    """Run *n_runs* independent simulations and return the total number of
    susceptible patients that became carriers across all runs.

    *make_carrier(run)* must return a fresh carrier Patient for each run.
    """
    total = 0
    for run in range(n_runs):
        cfg = config_factory()
        sim = MacroSimulator(config=cfg, n_hospitals=1, seed=run + seed_offset)
        sus = [
            Patient(patient_id=f"sus_{i}", state=HealthState.SUSCEPTIBLE)
            for i in range(n_susceptible)
        ]
        car = make_carrier(run)
        for p in sus + [car]:
            sim.admit(p, hospital_id=_H1, department=Department.WARD)
        for _ in range(n_days):
            sim.step()
        total += sum(1 for p in sus if p.state == HealthState.CARRIER)
    return total


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def susceptible_patient() -> Patient:
    """A patient in SUSCEPTIBLE state with default attributes."""
    return Patient(patient_id="sus_001", state=HealthState.SUSCEPTIBLE)


@pytest.fixture
def carrier_patient() -> Patient:
    """A patient in CARRIER state with all relevant fields populated."""
    return Patient(
        patient_id="car_001",
        state=HealthState.CARRIER,
        episode_id="episode_c001",
        age_years=60,
        compliance=0.85,
        vulnerability=1.2,
        immune_strength=0.9,
        immune_status="normal",
        sociability=1.0,
        resistant_fraction=0.3,
        dominant_genotype="R1",
        relative_transmissibility=1.5,
        p_clearance=0.05,
    )


@pytest.fixture
def isolated_carrier_patient() -> Patient:
    """A carrier who is currently isolated (context already set)."""
    patient = Patient(
        patient_id="iso_car_001",
        state=HealthState.CARRIER,
        episode_id="episode_iso001",
        age_years=70,
        compliance=0.9,
        vulnerability=1.5,
        immune_strength=0.7,
        immune_status="suppressed",
        sociability=0.5,
        resistant_fraction=0.6,
        dominant_genotype="R2",
        relative_transmissibility=1.8,
        p_clearance=0.03,
    )
    patient.update_context(
        PatientDailyContext(
            hospital_id=_H1,
            department=Department.WARD,
            hygiene_level=0.9,
            isolation_effectiveness=0.95,
            diagnostic_speed=0.8,
            is_isolated=True,
            regimen=AntibioticRegimen(on=True, abx_class="glycopeptide", dose_level="high"),
        )
    )
    return patient


@pytest.fixture
def icu_carrier_patient() -> Patient:
    """A carrier placed in the ICU (context already set)."""
    patient = Patient(
        patient_id="icu_car_001",
        state=HealthState.CARRIER,
        episode_id="episode_icu001",
        age_years=65,
        compliance=0.8,
        vulnerability=2.0,
        immune_strength=0.5,
        immune_status="suppressed",
        sociability=0.3,
        resistant_fraction=0.5,
        dominant_genotype="R1",
        relative_transmissibility=1.4,
        p_clearance=0.02,
    )
    patient.update_context(
        PatientDailyContext(
            hospital_id=_H1,
            department=Department.ICU,
            hygiene_level=0.95,
            isolation_effectiveness=0.9,
            diagnostic_speed=0.9,
            is_isolated=False,
            regimen=AntibioticRegimen(on=True, abx_class="beta_lactam", dose_level="high"),
        )
    )
    return patient


@pytest.fixture
def ward_carrier_patient() -> Patient:
    """A carrier placed in the general ward (context already set)."""
    patient = Patient(
        patient_id="ward_car_001",
        state=HealthState.CARRIER,
        episode_id="episode_ward001",
        age_years=55,
        compliance=0.7,
        vulnerability=1.0,
        immune_strength=1.0,
        immune_status="normal",
        sociability=1.0,
        resistant_fraction=0.2,
        dominant_genotype="R1",
        relative_transmissibility=1.0,
        p_clearance=0.05,
    )
    patient.update_context(
        PatientDailyContext(
            hospital_id=_H2,
            department=Department.WARD,
            hygiene_level=0.75,
            isolation_effectiveness=0.6,
            diagnostic_speed=0.4,
            is_isolated=False,
            regimen=AntibioticRegimen(on=False),
        )
    )
    return patient


@pytest.fixture
def default_config() -> SimulationConfig:
    """Default simulation configuration."""
    return SimulationConfig()


@pytest.fixture
def hospital(default_config: SimulationConfig) -> MacroSimulator:
    """Single-hospital macro simulator."""
    return MacroSimulator(config=default_config, n_hospitals=1, seed=42)


@pytest.fixture
def hospital_network(default_config: SimulationConfig) -> MacroSimulator:
    """10-hospital network macro simulator."""
    return MacroSimulator(config=default_config, n_hospitals=10, seed=42)


# =============================================================================
# Hospital / Placement Tests
# =============================================================================


class TestHospitalPlacement:
    """Patients can be placed in hospitals and departments."""

    def test_patient_assigned_to_hospital(
        self, hospital: MacroSimulator, susceptible_patient: Patient
    ):
        """Admitting assigns the patient to the specified hospital."""
        hospital.admit(susceptible_patient, hospital_id=_H1, department=Department.WARD)
        assert susceptible_patient.hospital_id == _H1

    @pytest.mark.parametrize("department", [Department.WARD, Department.ICU])
    def test_department_assignment(
        self,
        department: Department,
        hospital: MacroSimulator,
        susceptible_patient: Patient,
    ):
        """Patient is placed in the requested department."""
        hospital.admit(susceptible_patient, hospital_id=_H1, department=department)
        assert susceptible_patient.department == department

    def test_occupancy_increases_on_admit(
        self, hospital: MacroSimulator, susceptible_patient: Patient
    ):
        """Occupancy increments by one after admission."""
        before = hospital.get_occupancy(_H1)
        hospital.admit(susceptible_patient, hospital_id=_H1, department=Department.WARD)
        assert hospital.get_occupancy(_H1) == before + 1

    def test_occupancy_decreases_on_discharge(
        self, hospital: MacroSimulator, susceptible_patient: Patient
    ):
        """Occupancy decrements by one after discharge."""
        hospital.admit(susceptible_patient, hospital_id=_H1, department=Department.WARD)
        before = hospital.get_occupancy(_H1)
        hospital.discharge(susceptible_patient)
        assert hospital.get_occupancy(_H1) == before - 1


# =============================================================================
# Admission Tests
# =============================================================================


class TestAdmission:
    """The admission process works correctly."""

    def test_admit_adds_patient_to_registry(
        self, hospital: MacroSimulator, susceptible_patient: Patient
    ):
        """An admitted patient appears in the hospital's patient registry."""
        hospital.admit(susceptible_patient, hospital_id=_H1, department=Department.WARD)
        assert susceptible_patient in hospital.get_patients(_H1)

    @pytest.mark.parametrize("department", [Department.WARD, Department.ICU])
    def test_admit_to_department(
        self,
        department: Department,
        hospital: MacroSimulator,
        susceptible_patient: Patient,
    ):
        """Admission can target either ward or ICU."""
        hospital.admit(susceptible_patient, hospital_id=_H1, department=department)
        assert susceptible_patient.department == department

    def test_admitted_patient_gets_context_after_step(
        self, hospital: MacroSimulator, susceptible_patient: Patient
    ):
        """After one step the patient should have received a context update.

        We verify via public attributes that are written by update_context().
        """
        hospital.admit(susceptible_patient, hospital_id=_H1, department=Department.WARD)
        hospital.step()
        # Public attributes set by update_context
        assert susceptible_patient.hospital_id == _H1
        assert susceptible_patient.department == Department.WARD
        assert isinstance(susceptible_patient.regimen, AntibioticRegimen)

    def test_admit_susceptible_no_carrier_artifacts(
        self, hospital: MacroSimulator, susceptible_patient: Patient
    ):
        """Admitting a susceptible patient must not create carrier-specific state."""
        hospital.admit(susceptible_patient, hospital_id=_H1, department=Department.WARD)
        assert susceptible_patient.state == HealthState.SUSCEPTIBLE
        assert susceptible_patient.episode_id is None
        assert susceptible_patient.resistant_fraction == 0.0
        assert susceptible_patient.dominant_genotype == "S"

    def test_admit_carrier_preserves_episode(
        self, hospital: MacroSimulator, carrier_patient: Patient
    ):
        """Admitting a carrier preserves episode_id and carrier state."""
        original_episode = carrier_patient.episode_id
        hospital.admit(carrier_patient, hospital_id=_H1, department=Department.WARD)
        assert carrier_patient.state == HealthState.CARRIER
        assert carrier_patient.episode_id == original_episode


# =============================================================================
# Discharge Tests
# =============================================================================


class TestDischarge:
    """The discharge process works correctly."""

    def test_discharge_removes_from_occupancy(
        self, hospital: MacroSimulator, susceptible_patient: Patient
    ):
        """Discharging drops occupancy to zero (single-patient hospital)."""
        hospital.admit(susceptible_patient, hospital_id=_H1, department=Department.WARD)
        hospital.discharge(susceptible_patient)
        assert hospital.get_occupancy(_H1) == 0

    def test_discharged_not_in_registry(
        self, hospital: MacroSimulator, susceptible_patient: Patient
    ):
        """A discharged patient no longer appears in the registry."""
        hospital.admit(susceptible_patient, hospital_id=_H1, department=Department.WARD)
        hospital.discharge(susceptible_patient)
        assert susceptible_patient not in hospital.get_patients(_H1)

    def test_discharge_clears_hospital_id(self, hospital: MacroSimulator, carrier_patient: Patient):
        """Discharge preserves the Patient object but sets hospital_id to None."""
        hospital.admit(carrier_patient, hospital_id=_H1, department=Department.WARD)
        pid = carrier_patient.patient_id
        hospital.discharge(carrier_patient)
        assert carrier_patient.patient_id == pid
        assert carrier_patient.hospital_id is None

    def test_carrier_discharge_preserves_episode_state(
        self, hospital: MacroSimulator, carrier_patient: Patient
    ):
        """Discharging a carrier must not corrupt episode / micro-derived state."""
        snap = {
            "episode_id": carrier_patient.episode_id,
            "resistant_fraction": carrier_patient.resistant_fraction,
            "dominant_genotype": carrier_patient.dominant_genotype,
        }
        hospital.admit(carrier_patient, hospital_id=_H1, department=Department.WARD)
        hospital.discharge(carrier_patient)
        assert carrier_patient.state == HealthState.CARRIER
        assert carrier_patient.episode_id == snap["episode_id"]
        assert carrier_patient.resistant_fraction == snap["resistant_fraction"]
        assert carrier_patient.dominant_genotype == snap["dominant_genotype"]

    def test_discharge_non_admitted_raises(
        self, hospital: MacroSimulator, susceptible_patient: Patient
    ):
        """Discharging a patient who is not admitted should raise ValueError."""
        with pytest.raises(ValueError):
            hospital.discharge(susceptible_patient)


# =============================================================================
# Transfer Tests
# =============================================================================


class TestTransfer:
    """Inter-hospital transfers work correctly."""

    def test_transfer_moves_patient(
        self, hospital_network: MacroSimulator, susceptible_patient: Patient
    ):
        """Transfer moves the patient to the destination hospital."""
        hospital_network.admit(susceptible_patient, hospital_id=_H1, department=Department.WARD)
        hospital_network.transfer(susceptible_patient, to_hospital_id=_H2)
        assert susceptible_patient.hospital_id == _H2

    def test_transfer_preserves_all_carrier_fields(
        self, hospital_network: MacroSimulator, carrier_patient: Patient
    ):
        """Transfer must not mutate any patient identity or micro-derived field."""
        snap = {
            "patient_id": carrier_patient.patient_id,
            "state": carrier_patient.state,
            "episode_id": carrier_patient.episode_id,
            "resistant_fraction": carrier_patient.resistant_fraction,
            "dominant_genotype": carrier_patient.dominant_genotype,
            "relative_transmissibility": carrier_patient.relative_transmissibility,
            "p_clearance": carrier_patient.p_clearance,
        }
        hospital_network.admit(carrier_patient, hospital_id=_H1, department=Department.WARD)
        hospital_network.transfer(carrier_patient, to_hospital_id=_H2)

        for field, expected in snap.items():
            assert getattr(carrier_patient, field) == expected, (
                f"Transfer mutated {field}: expected {expected}, "
                f"got {getattr(carrier_patient, field)}"
            )

    def test_transfer_updates_both_occupancies(
        self, hospital_network: MacroSimulator, susceptible_patient: Patient
    ):
        """Source occupancy decrements, destination increments."""
        hospital_network.admit(susceptible_patient, hospital_id=_H1, department=Department.WARD)
        src_before = hospital_network.get_occupancy(_H1)
        dst_before = hospital_network.get_occupancy(_H2)
        hospital_network.transfer(susceptible_patient, to_hospital_id=_H2)
        assert hospital_network.get_occupancy(_H1) == src_before - 1
        assert hospital_network.get_occupancy(_H2) == dst_before + 1

    @pytest.mark.parametrize(
        "state, episode_id",
        [
            (HealthState.SUSCEPTIBLE, None),
            (HealthState.CARRIER, "ep_xfer"),
        ],
        ids=["susceptible", "carrier"],
    )
    def test_transfer_works_for_both_states(
        self,
        state: HealthState,
        episode_id: str | None,
        hospital_network: MacroSimulator,
    ):
        """Both susceptible and carrier patients can be transferred."""
        patient = Patient(patient_id="xfer_p", state=state, episode_id=episode_id)
        hospital_network.admit(patient, hospital_id=_H1, department=Department.WARD)
        hospital_network.transfer(patient, to_hospital_id=_H3)
        assert patient.state == state
        assert patient.hospital_id == _H3

    def test_transfer_chain_across_network(
        self, hospital_network: MacroSimulator, carrier_patient: Patient
    ):
        """Repeated transfers across all 10 hospitals preserve state."""
        hospital_network.admit(carrier_patient, hospital_id=_H1, department=Department.WARD)
        for i in range(2, 11):
            hid = f"hospital_{i:03d}"
            hospital_network.transfer(carrier_patient, to_hospital_id=hid)
            hospital_network.step()
            assert carrier_patient.hospital_id == hid
            assert carrier_patient.state == HealthState.CARRIER
            assert carrier_patient.episode_id == "episode_c001"

    def test_transfer_non_admitted_raises(
        self, hospital_network: MacroSimulator, susceptible_patient: Patient
    ):
        """Transferring a non-admitted patient should raise ValueError."""
        with pytest.raises(ValueError):
            hospital_network.transfer(susceptible_patient, to_hospital_id=_H2)

    def test_transfer_to_same_hospital_is_safe(
        self, hospital_network: MacroSimulator, susceptible_patient: Patient
    ):
        """Transfer to the same hospital is either a no-op or raises ValueError."""
        hospital_network.admit(susceptible_patient, hospital_id=_H1, department=Department.WARD)
        occ_before = hospital_network.get_occupancy(_H1)
        try:
            hospital_network.transfer(susceptible_patient, to_hospital_id=_H1)
        except ValueError:
            pass  # explicitly forbidden — acceptable
        else:
            # no-op — occupancy and placement unchanged
            assert susceptible_patient.hospital_id == _H1
            assert hospital_network.get_occupancy(_H1) == occ_before


# =============================================================================
# Daily Context Update Tests
# =============================================================================


class TestDailyContextUpdate:
    """The macro layer updates PatientDailyContext correctly each day."""

    def test_public_fields_set_after_step(
        self, hospital: MacroSimulator, susceptible_patient: Patient
    ):
        """After a step, the public fields written by update_context() must be set."""
        hospital.admit(susceptible_patient, hospital_id=_H1, department=Department.WARD)
        hospital.step()
        assert susceptible_patient.hospital_id == _H1
        assert susceptible_patient.department == Department.WARD
        assert isinstance(susceptible_patient.is_isolated, bool)
        assert isinstance(susceptible_patient.regimen, AntibioticRegimen)

    def test_context_contains_environment_fields(
        self, hospital: MacroSimulator, susceptible_patient: Patient
    ):
        """The stored PatientDailyContext must carry valid environment fields.

        hygiene_level, isolation_effectiveness, and diagnostic_speed live only
        on PatientDailyContext (not mirrored on Patient), so we access _ctx here.
        """
        hospital.admit(susceptible_patient, hospital_id=_H1, department=Department.WARD)
        hospital.step()
        ctx = susceptible_patient._ctx
        assert ctx is not None
        assert 0.0 <= ctx.hygiene_level <= 1.0
        assert 0.0 <= ctx.isolation_effectiveness <= 1.0
        assert ctx.diagnostic_speed >= 0.0

    @pytest.mark.parametrize("department", [Department.WARD, Department.ICU])
    def test_context_department_matches_placement(
        self, department: Department, hospital: MacroSimulator
    ):
        """Context department must match the department the patient was placed in."""
        patient = Patient(patient_id=f"dept_{department.value}", state=HealthState.SUSCEPTIBLE)
        hospital.admit(patient, hospital_id=_H1, department=department)
        hospital.step()
        assert patient.department == department


# =============================================================================
# Transmission Tests
# =============================================================================


class TestTransmission:
    """Within-hospital transmission behaviour.

    Transmission is stochastic, so we test invariants and monotonic relationships
    over many Monte-Carlo repetitions rather than exact outcomes.
    """

    MC_RUNS = 50
    MC_DAYS = 10

    def test_only_susceptible_can_become_carrier(self, hospital: MacroSimulator):
        """Only patients starting as SUSCEPTIBLE may transition to CARRIER."""
        patients = _make_patients(n_susceptible=5, n_carrier=3)
        for p in patients:
            hospital.admit(p, hospital_id=_H1, department=Department.WARD)
        for _ in range(30):
            hospital.step()
        for p in patients:
            if p.state == HealthState.CARRIER and p.patient_id.startswith("sus_"):
                assert p.episode_id is not None

    def test_more_carriers_more_colonization(self):
        """More seed carriers ⇒ at least as many new colonizations (monotonic)."""
        low = _count_new_carriers_mc(
            n_runs=self.MC_RUNS,
            n_days=self.MC_DAYS,
            n_susceptible=10,
            make_carrier=lambda run: Patient(
                patient_id="car_0",
                state=HealthState.CARRIER,
                episode_id="ep0",
                p_clearance=0.0,
            ),
        )
        # For "high" we place 5 carriers; easiest via a direct loop
        high = 0
        for run in range(self.MC_RUNS):
            sim = MacroSimulator(config=SimulationConfig(), n_hospitals=1, seed=run + 10_000)
            patients = _make_patients(n_susceptible=10, n_carrier=5, carrier_p_clearance=0.0)
            for p in patients:
                sim.admit(p, hospital_id=_H1, department=Department.WARD)
            for _ in range(self.MC_DAYS):
                sim.step()
            high += sum(
                1
                for p in patients
                if p.state == HealthState.CARRIER and p.patient_id.startswith("sus_")
            )
        assert high >= low, f"More carriers should produce >= colonizations: high={high}, low={low}"

    def test_higher_transmissibility_higher_pressure(self):
        """Higher relative_transmissibility ⇒ at least as many colonizations."""
        low = _count_new_carriers_mc(
            n_runs=self.MC_RUNS,
            n_days=self.MC_DAYS,
            n_susceptible=10,
            make_carrier=lambda _: Patient(
                patient_id="car_lt",
                state=HealthState.CARRIER,
                episode_id="ep_lt",
                relative_transmissibility=0.5,
                p_clearance=0.0,
            ),
        )
        high = _count_new_carriers_mc(
            n_runs=self.MC_RUNS,
            n_days=self.MC_DAYS,
            n_susceptible=10,
            make_carrier=lambda _: Patient(
                patient_id="car_ht",
                state=HealthState.CARRIER,
                episode_id="ep_ht",
                relative_transmissibility=3.0,
                p_clearance=0.0,
            ),
        )
        assert high >= low

    def test_isolation_reduces_transmission(self):
        """An isolated carrier should cause fewer colonizations than a non-isolated one."""
        no_iso = _count_new_carriers_mc(
            n_runs=self.MC_RUNS,
            n_days=self.MC_DAYS,
            n_susceptible=10,
            make_carrier=lambda _: Patient(
                patient_id="car",
                state=HealthState.CARRIER,
                episode_id="ep",
                relative_transmissibility=2.0,
                p_clearance=0.0,
                is_isolated=False,
            ),
        )
        iso = _count_new_carriers_mc(
            n_runs=self.MC_RUNS,
            n_days=self.MC_DAYS,
            n_susceptible=10,
            make_carrier=lambda _: Patient(
                patient_id="car_iso",
                state=HealthState.CARRIER,
                episode_id="ep_iso",
                relative_transmissibility=2.0,
                p_clearance=0.0,
                is_isolated=True,
            ),
        )
        assert iso <= no_iso, f"Isolation should reduce colonizations: iso={iso}, no_iso={no_iso}"

    def test_higher_hygiene_reduces_transmission(self):
        """Higher base_hygiene ⇒ fewer colonizations."""
        low_hyg = _count_new_carriers_mc(
            n_runs=self.MC_RUNS,
            n_days=self.MC_DAYS,
            n_susceptible=10,
            make_carrier=lambda _: Patient(
                patient_id="car",
                state=HealthState.CARRIER,
                episode_id="ep",
                p_clearance=0.0,
            ),
            config_factory=lambda: SimulationConfig(base_hygiene=0.2),
        )
        high_hyg = _count_new_carriers_mc(
            n_runs=self.MC_RUNS,
            n_days=self.MC_DAYS,
            n_susceptible=10,
            make_carrier=lambda _: Patient(
                patient_id="car",
                state=HealthState.CARRIER,
                episode_id="ep",
                p_clearance=0.0,
            ),
            config_factory=lambda: SimulationConfig(base_hygiene=0.95),
            seed_offset=0,
        )
        assert high_hyg <= low_hyg

    def test_no_carriers_no_colonizations(self, hospital: MacroSimulator):
        """Without any carriers no susceptible patient can become colonised."""
        patients = [
            Patient(patient_id=f"sus_{i}", state=HealthState.SUSCEPTIBLE) for i in range(20)
        ]
        for p in patients:
            hospital.admit(p, hospital_id=_H1, department=Department.WARD)
        for _ in range(50):
            hospital.step()
        for p in patients:
            assert (
                p.state == HealthState.SUSCEPTIBLE
            ), f"{p.patient_id} became carrier with no source carriers"

    def test_only_carriers_no_new_colonizations(self, hospital: MacroSimulator):
        """If every patient is already a carrier (clearance=0), no new colonizations occur."""
        patients = [
            Patient(
                patient_id=f"car_{i}",
                state=HealthState.CARRIER,
                episode_id=f"ep_{i}",
                p_clearance=0.0,
            )
            for i in range(10)
        ]
        for p in patients:
            hospital.admit(p, hospital_id=_H1, department=Department.WARD)
        for _ in range(20):
            hospital.step()
        for p in patients:
            assert p.state == HealthState.CARRIER


# =============================================================================
# Clearance Tests
# =============================================================================


class TestClearance:
    """Macro handles C→S clearance correctly."""

    def test_susceptible_never_clears(self, hospital: MacroSimulator, susceptible_patient: Patient):
        """Clearance must be a no-op for susceptible patients."""
        hospital.admit(susceptible_patient, hospital_id=_H1, department=Department.WARD)
        for _ in range(50):
            hospital.step()
        assert susceptible_patient.state == HealthState.SUSCEPTIBLE

    def test_p_clearance_1_always_clears(self, hospital: MacroSimulator):
        """A carrier with p_clearance=1.0 must always clear in one day."""
        patient = Patient(
            patient_id="clear_1",
            state=HealthState.CARRIER,
            episode_id="ep_c1",
            p_clearance=1.0,
        )
        hospital.admit(patient, hospital_id=_H1, department=Department.WARD)
        hospital.step()
        assert patient.state == HealthState.SUSCEPTIBLE

    def test_p_clearance_0_never_clears(self, hospital: MacroSimulator):
        """A carrier with p_clearance=0.0 must never clear."""
        patient = Patient(
            patient_id="clear_0",
            state=HealthState.CARRIER,
            episode_id="ep_c0",
            p_clearance=0.0,
        )
        hospital.admit(patient, hospital_id=_H1, department=Department.WARD)
        hospital.step()
        assert patient.state == HealthState.CARRIER

    def test_clearance_resets_all_carrier_fields(self, hospital: MacroSimulator):
        """After clearance, all carrier-specific fields must match clear_carriage() defaults."""
        patient = Patient(
            patient_id="clear_reset",
            state=HealthState.CARRIER,
            episode_id="ep_reset",
            resistant_fraction=0.8,
            dominant_genotype="R2",
            relative_transmissibility=2.0,
            p_clearance=1.0,
            severity_modifier=1.5,
            lethality_modifier=1.3,
        )
        hospital.admit(patient, hospital_id=_H1, department=Department.WARD)
        hospital.step()

        assert patient.state == HealthState.SUSCEPTIBLE
        assert patient.episode_id is None
        assert patient.resistant_fraction == 0.0
        assert patient.dominant_genotype == "S"
        assert patient.relative_transmissibility == 1.0
        assert patient.p_clearance == 0.02
        assert patient.severity_modifier == 1.0
        assert patient.lethality_modifier == 1.0

    def test_p_clearance_respected_statistically(self):
        """With p_clearance=0.5, roughly half of 200 runs should clear after one day."""
        n_runs = 200
        cleared = 0
        for run in range(n_runs):
            sim = MacroSimulator(config=SimulationConfig(), n_hospitals=1, seed=run)
            patient = Patient(
                patient_id="clr_half",
                state=HealthState.CARRIER,
                episode_id=f"ep_{run}",
                p_clearance=0.5,
            )
            sim.admit(patient, hospital_id=_H1, department=Department.WARD)
            sim.step()
            if patient.state == HealthState.SUSCEPTIBLE:
                cleared += 1
        # Binomial(200, 0.5): 99.9 % CI ≈ [72, 128]
        assert 60 < cleared < 140, f"Expected ~100 clearances, got {cleared}"


# =============================================================================
# Antibiotic Policy Tests
# =============================================================================


class TestAntibioticPolicy:
    """Department-specific antibiotic policies."""

    def test_icu_higher_abx_frequency_than_ward(self):
        """ICU patients should receive antibiotics more often than ward patients.

        Measured as fraction of patient-days with regimen.on=True over 100 runs.
        Uses the public patient.regimen attribute (set by update_context).
        """
        n_runs = 100
        n_days = 20
        icu_on = ward_on = 0
        total_days = n_runs * n_days  # same denominator for both

        for run in range(n_runs):
            sim = MacroSimulator(config=SimulationConfig(), n_hospitals=1, seed=run)
            icu_p = Patient(patient_id=f"icu_{run}", state=HealthState.SUSCEPTIBLE)
            ward_p = Patient(patient_id=f"ward_{run}", state=HealthState.SUSCEPTIBLE)
            sim.admit(icu_p, hospital_id=_H1, department=Department.ICU)
            sim.admit(ward_p, hospital_id=_H1, department=Department.WARD)
            for _ in range(n_days):
                sim.step()
                if icu_p.regimen.on:
                    icu_on += 1
                if ward_p.regimen.on:
                    ward_on += 1

        icu_frac = icu_on / total_days
        ward_frac = ward_on / total_days
        assert (
            icu_frac > ward_frac
        ), f"ICU ABX fraction ({icu_frac:.3f}) should exceed ward ({ward_frac:.3f})"

    def test_regimen_types_valid_after_step(
        self, hospital: MacroSimulator, susceptible_patient: Patient
    ):
        """The regimen set by macro must be a well-typed AntibioticRegimen."""
        hospital.admit(susceptible_patient, hospital_id=_H1, department=Department.WARD)
        hospital.step()
        reg = susceptible_patient.regimen
        assert isinstance(reg, AntibioticRegimen)
        assert isinstance(reg.on, bool)
        assert isinstance(reg.abx_class, str)
        assert isinstance(reg.dose_level, str)

    def test_no_abx_implies_class_none(self, hospital: MacroSimulator):
        """When regimen.on is False the abx_class should be 'none'."""
        patient = Patient(patient_id="no_abx", state=HealthState.SUSCEPTIBLE)
        hospital.admit(patient, hospital_id=_H1, department=Department.WARD)
        hospital.step()
        reg = patient.regimen
        if not reg.on:
            assert reg.abx_class == "none"

    def test_abx_on_implies_valid_class_and_dose(self, hospital: MacroSimulator):
        """When regimen.on is True, abx_class != 'none' and dose_level is valid."""
        patient = Patient(
            patient_id="abx_on",
            state=HealthState.CARRIER,
            episode_id="ep_abx",
        )
        hospital.admit(patient, hospital_id=_H1, department=Department.ICU)
        for _ in range(30):
            hospital.step()
            reg = patient.regimen
            if reg.on:
                assert reg.abx_class != "none" and len(reg.abx_class) > 0
                assert reg.dose_level in {"low", "std", "high"}
                return
        pytest.skip("No ABX-on day observed in 30 steps for ICU carrier")


# =============================================================================
# Isolation and Detection Tests
# =============================================================================


class TestIsolationAndDetection:
    """Isolation and diagnostic detection are handled correctly."""

    def test_isolated_flag_is_bool_after_step(self, hospital: MacroSimulator):
        """After a step the isolation flag must be a boolean."""
        patient = Patient(
            patient_id="det",
            state=HealthState.CARRIER,
            episode_id="ep_det",
            is_isolated=True,
        )
        hospital.admit(patient, hospital_id=_H1, department=Department.WARD)
        hospital.step()
        assert isinstance(patient.is_isolated, bool)

    def test_diagnostic_speed_available(
        self, hospital: MacroSimulator, susceptible_patient: Patient
    ):
        """diagnostic_speed must be a non-negative float in the daily context."""
        hospital.admit(susceptible_patient, hospital_id=_H1, department=Department.WARD)
        hospital.step()
        ctx = susceptible_patient._ctx
        assert isinstance(ctx.diagnostic_speed, (int, float))
        assert ctx.diagnostic_speed >= 0.0


# =============================================================================
# Multi-Day Simulation Tests
# =============================================================================


class TestMultiDaySimulation:
    """Multi-day simulations are stable and consistent."""

    def test_carrier_persists_across_days(self, hospital: MacroSimulator, carrier_patient: Patient):
        """A carrier (p_clearance=0) must persist without state corruption."""
        carrier_patient.p_clearance = 0.0
        hospital.admit(carrier_patient, hospital_id=_H1, department=Department.WARD)
        for _ in range(10):
            hospital.step()
            assert carrier_patient.patient_id == "car_001"
            assert carrier_patient.state == HealthState.CARRIER
            assert carrier_patient.episode_id == "episode_c001"

    def test_hospital_assignment_stable(
        self, hospital: MacroSimulator, susceptible_patient: Patient
    ):
        """Hospital/department assignment stays consistent across days."""
        hospital.admit(susceptible_patient, hospital_id=_H1, department=Department.WARD)
        for _ in range(10):
            hospital.step()
            assert susceptible_patient.hospital_id == _H1
            assert susceptible_patient.department == Department.WARD

    def test_context_refreshed_every_day(
        self, hospital: MacroSimulator, susceptible_patient: Patient
    ):
        """update_context is called each step (regimen is always set)."""
        hospital.admit(susceptible_patient, hospital_id=_H1, department=Department.WARD)
        for _ in range(5):
            hospital.step()
            assert isinstance(susceptible_patient.regimen, AntibioticRegimen)

    def test_multiday_multi_hospital_no_crash(self, hospital_network: MacroSimulator):
        """30 days across 10 hospitals with mixed patients must not crash."""
        patients = _make_patients(n_susceptible=30, n_carrier=10, carrier_p_clearance=0.0)
        for i, p in enumerate(patients):
            hid = _HOSPITAL_IDS[i % 10]
            dept = Department.ICU if i % 5 == 0 else Department.WARD
            hospital_network.admit(p, hospital_id=hid, department=dept)
        for _ in range(30):
            hospital_network.step()
        total = sum(hospital_network.get_occupancy(h) for h in _HOSPITAL_IDS)
        assert total > 0


# =============================================================================
# Batch / Network Tests
# =============================================================================


class TestBatchAndNetwork:
    """Batch processing of many patients across a hospital network."""

    def test_occupancy_across_10_hospitals(self, hospital_network: MacroSimulator):
        """Total occupancy equals the number of admitted patients after one step."""
        patients = _make_patients(n_susceptible=50, n_carrier=20)
        for i, p in enumerate(patients):
            hospital_network.admit(p, hospital_id=_HOSPITAL_IDS[i % 10], department=Department.WARD)
        hospital_network.step()
        total = sum(hospital_network.get_occupancy(h) for h in _HOSPITAL_IDS)
        assert total == len(patients)

    def test_mixed_population_valid_states(self, hospital_network: MacroSimulator):
        """After simulation, every patient must still be in a valid macro state."""
        patients = _make_patients(n_susceptible=15, n_carrier=5)
        for p in patients:
            hospital_network.admit(p, hospital_id=_H1, department=Department.WARD)
        for _ in range(10):
            hospital_network.step()
        for p in patients:
            assert p.state in _VALID_STATES

    def test_random_transfers_over_time(self, hospital_network: MacroSimulator):
        """Random daily transfers across the network must not corrupt state."""
        patients = _make_patients(n_susceptible=10, n_carrier=5)
        for i, p in enumerate(patients):
            hospital_network.admit(p, hospital_id=_HOSPITAL_IDS[i % 10], department=Department.WARD)
        rng = random.Random(99)
        for _ in range(10):
            admitted = [p for p in patients if p.hospital_id is not None]
            if admitted:
                p = rng.choice(admitted)
                dest = rng.choice(_HOSPITAL_IDS)
                if dest != p.hospital_id:
                    hospital_network.transfer(p, to_hospital_id=dest)
            hospital_network.step()
        for p in patients:
            assert p.state in _VALID_STATES

    def test_patient_count_conserved(self, hospital_network: MacroSimulator):
        """Without admissions/discharges, total count is conserved."""
        patients = _make_patients(n_susceptible=30, n_carrier=10)
        for i, p in enumerate(patients):
            hospital_network.admit(p, hospital_id=_HOSPITAL_IDS[i % 10], department=Department.WARD)
        for _ in range(20):
            hospital_network.step()
        total = sum(hospital_network.get_occupancy(h) for h in _HOSPITAL_IDS)
        assert total == len(patients)

    def test_patient_ids_stable(self, hospital_network: MacroSimulator):
        """Patient IDs must remain unique and unchanged."""
        patients = _make_patients(n_susceptible=20, n_carrier=10)
        original_ids = {p.patient_id for p in patients}
        for p in patients:
            hospital_network.admit(p, hospital_id=_H1, department=Department.WARD)
        for _ in range(10):
            hospital_network.step()
        assert {p.patient_id for p in patients} == original_ids


# =============================================================================
# Edge Case Tests
# =============================================================================


class TestEdgeCases:
    """Edge cases and boundary conditions."""

    def test_empty_hospital_step(self, hospital: MacroSimulator):
        """Stepping an empty hospital must not crash."""
        hospital.step()
        assert hospital.get_occupancy(_H1) == 0

    def test_zero_sociability_no_spread(self, hospital: MacroSimulator):
        """A carrier with sociability=0 contributes zero transmission (multiplier is 0)."""
        sus = [Patient(patient_id=f"s_{i}", state=HealthState.SUSCEPTIBLE) for i in range(10)]
        car = Patient(
            patient_id="zero_soc",
            state=HealthState.CARRIER,
            episode_id="ep_zs",
            sociability=0.0,
            p_clearance=0.0,
        )
        for p in sus + [car]:
            hospital.admit(p, hospital_id=_H1, department=Department.WARD)
        for _ in range(30):
            hospital.step()
        for p in sus:
            assert (
                p.state == HealthState.SUSCEPTIBLE
            ), f"{p.patient_id} colonized despite zero-sociability carrier"

    @pytest.mark.parametrize(
        "label, kwargs",
        [
            ("high_vulnerability", dict(vulnerability=3.0, immune_strength=0.2)),
            (
                "immunosuppressed",
                dict(
                    immune_status="suppressed",
                    immune_strength=0.3,
                    vulnerability=2.0,
                    history_flags={"immunosuppression"},
                ),
            ),
            ("zero_compliance", dict(compliance=0.0)),
        ],
        ids=["high_vulnerability", "immunosuppressed", "zero_compliance"],
    )
    def test_extreme_patient_processable(self, label: str, kwargs: dict, hospital: MacroSimulator):
        """Patients with extreme attribute values must be processable."""
        patient = Patient(patient_id=label, state=HealthState.SUSCEPTIBLE, **kwargs)
        hospital.admit(patient, hospital_id=_H1, department=Department.ICU)
        for _ in range(10):
            hospital.step()
        assert patient.state in _VALID_STATES

    def test_missing_context_filled_by_macro(self, hospital: MacroSimulator):
        """A patient with no prior context receives one from macro after step."""
        patient = Patient(patient_id="no_ctx", state=HealthState.SUSCEPTIBLE)
        assert patient._ctx is None
        hospital.admit(patient, hospital_id=_H1, department=Department.WARD)
        hospital.step()
        # After step, public fields prove context was set
        assert patient.hospital_id == _H1
        assert isinstance(patient.regimen, AntibioticRegimen)

    def test_isolated_carrier_minimal_spread(self, hospital: MacroSimulator):
        """An isolated carrier should not crash the simulation."""
        sus = [Patient(patient_id=f"s_{i}", state=HealthState.SUSCEPTIBLE) for i in range(10)]
        car = Patient(
            patient_id="iso_car",
            state=HealthState.CARRIER,
            episode_id="ep_iso",
            p_clearance=0.0,
            is_isolated=True,
        )
        for p in sus + [car]:
            hospital.admit(p, hospital_id=_H1, department=Department.WARD)
        for _ in range(10):
            hospital.step()
        for p in sus:
            assert p.state in _VALID_STATES


# =============================================================================
# Reproducibility Tests
# =============================================================================


class TestReproducibility:
    """Macro simulations are reproducible with the same seed."""

    @staticmethod
    def _run_sim(seed: int) -> list[str]:
        sim = MacroSimulator(config=SimulationConfig(), n_hospitals=1, seed=seed)
        patients = _make_patients(n_susceptible=10, n_carrier=3)
        for p in patients:
            sim.admit(p, hospital_id=_H1, department=Department.WARD)
        for _ in range(20):
            sim.step()
        return [p.state.value for p in patients]

    def test_same_seed_same_outcome(self):
        """Identical seeds must produce identical state sequences."""
        assert self._run_sim(12345) == self._run_sim(12345)

    def test_different_seeds_diverge(self):
        """20 different seeds should yield at least two distinct outcomes."""
        results = {tuple(self._run_sim(s)) for s in range(20)}
        assert len(results) > 1, "All 20 seeds produced identical outcomes"


# =============================================================================
# Patient Interface Compatibility Tests
# =============================================================================


class TestPatientInterfaceCompatibility:
    """Macro works entirely through the Patient public interface and does not
    break micro-layer assumptions.
    """

    def test_macro_runs_without_micro(self, hospital: MacroSimulator, carrier_patient: Patient):
        """Macro completes without importing or running the micro simulator."""
        hospital.admit(carrier_patient, hospital_id=_H1, department=Department.WARD)
        for _ in range(5):
            hospital.step()

    def test_macro_step_calls_micro_once_per_carrier_per_day(self, hospital: MacroSimulator):
        """When micro is provided, each active carrier is processed once per macro day."""

        class _RecordingMicro:
            def __init__(self):
                self.calls = 0
                self.request_batches: list[list[dict]] = []

            def process_batch(self, requests: list[dict], parallel: bool = True) -> list[dict]:
                self.calls += 1
                self.request_batches.append(list(requests))
                return [
                    {
                        "episode_id": req["episode_id"],
                        "patient_id": req["patient_id"],
                        "updated_state": {
                            "resistant_fraction": req["initial_state"]["resistant_fraction"],
                            "dominant_genotype": req["initial_state"]["dominant_genotype"],
                        },
                        "derived_effects": {
                            "relative_transmissibility": 1.0,
                            "p_clearance": 0.0,
                            "severity_modifier": 1.0,
                            "lethality_modifier": 1.0,
                        },
                    }
                    for req in requests
                ]

        micro = _RecordingMicro()
        carrier_a = Patient(patient_id="car_a", state=HealthState.CARRIER, episode_id="ep_a")
        carrier_b = Patient(patient_id="car_b", state=HealthState.CARRIER, episode_id="ep_b")
        susceptible = Patient(patient_id="sus_a", state=HealthState.SUSCEPTIBLE)

        for p in (carrier_a, carrier_b, susceptible):
            hospital.admit(p, hospital_id=_H1, department=Department.WARD)

        hospital.step(micro_simulator=micro, run_id="run_micro_phase")

        assert micro.calls == 1
        assert len(micro.request_batches[0]) == 2
        assert {req["patient_id"] for req in micro.request_batches[0]} == {"car_a", "car_b"}
        assert all(req["dt_days"] == 1 for req in micro.request_batches[0])

    def test_macro_batches_carriers_across_hospitals_once_per_day(self):
        """All active carriers for a day should be sent in one global micro batch."""

        class _RecordingMicro:
            def __init__(self):
                self.calls = 0
                self.batch_sizes: list[int] = []

            def process_batch(self, requests: list[dict], parallel: bool = True) -> list[dict]:
                self.calls += 1
                self.batch_sizes.append(len(requests))
                return [
                    {
                        "episode_id": req["episode_id"],
                        "patient_id": req["patient_id"],
                        "updated_state": {
                            "resistant_fraction": req["initial_state"]["resistant_fraction"],
                            "dominant_genotype": req["initial_state"]["dominant_genotype"],
                        },
                        "derived_effects": {
                            "relative_transmissibility": 1.0,
                            "p_clearance": 0.0,
                            "severity_modifier": 1.0,
                            "lethality_modifier": 1.0,
                        },
                    }
                    for req in requests
                ]

        sim = MacroSimulator(config=SimulationConfig(), n_hospitals=2, seed=42)
        sim.admit(
            Patient(patient_id="car_h1", state=HealthState.CARRIER, episode_id="ep_h1"),
            hospital_id=_H1,
            department=Department.WARD,
        )
        sim.admit(
            Patient(patient_id="car_h2", state=HealthState.CARRIER, episode_id="ep_h2"),
            hospital_id=_H2,
            department=Department.WARD,
        )

        micro = _RecordingMicro()
        sim.step(micro_simulator=micro, run_id="run_multi_hospital_batch")

        assert micro.calls == 1
        assert micro.batch_sizes == [2]

    def test_macro_applies_micro_response_before_next_clearance(self, hospital: MacroSimulator):
        """Micro-updated p_clearance should affect the next macro-day clearance phase."""

        class _ClearanceMicro:
            def process_batch(self, requests: list[dict], parallel: bool = True) -> list[dict]:
                return [
                    {
                        "episode_id": req["episode_id"],
                        "patient_id": req["patient_id"],
                        "updated_state": {
                            "resistant_fraction": req["initial_state"]["resistant_fraction"],
                            "dominant_genotype": req["initial_state"]["dominant_genotype"],
                        },
                        "derived_effects": {
                            "relative_transmissibility": 1.0,
                            "p_clearance": 1.0,
                            "severity_modifier": 1.0,
                            "lethality_modifier": 1.0,
                        },
                    }
                    for req in requests
                ]

        patient = Patient(
            patient_id="carrier_clear",
            state=HealthState.CARRIER,
            episode_id="ep_clear",
            p_clearance=0.0,
        )
        hospital.admit(patient, hospital_id=_H1, department=Department.WARD)

        micro = _ClearanceMicro()
        hospital.step(micro_simulator=micro, run_id="run_clear")
        assert patient.state == HealthState.CARRIER
        assert patient.p_clearance == 1.0

        hospital.step(micro_simulator=micro, run_id="run_clear")
        assert patient.state == HealthState.SUSCEPTIBLE

    def test_macro_clearance_removes_micro_episode_state(self, hospital: MacroSimulator):
        """Macro-owned clearance must also drop the persisted micro episode."""

        class _TrackingMicro:
            def __init__(self):
                self.cleared_episode_ids: list[str] = []

            def clear_episode(self, episode_id: str) -> None:
                self.cleared_episode_ids.append(episode_id)

            def process_batch(self, requests: list[dict], parallel: bool = True) -> list[dict]:
                raise AssertionError("Cleared episodes must not be sent to micro.")

        patient = Patient(
            patient_id="carrier_clear_cleanup",
            state=HealthState.CARRIER,
            episode_id="ep_cleanup",
            p_clearance=1.0,
        )
        hospital.admit(patient, hospital_id=_H1, department=Department.WARD)

        micro = _TrackingMicro()
        hospital.step(micro_simulator=micro, run_id="run_cleanup")

        assert patient.state == HealthState.SUSCEPTIBLE
        assert micro.cleared_episode_ids == ["ep_cleanup"]

    def test_micro_request_valid_after_macro_days(self, hospital: MacroSimulator):
        """After several macro days a carrier still produces a valid micro request."""
        patient = Patient(
            patient_id="compat_001",
            state=HealthState.CARRIER,
            episode_id="ep_compat",
            p_clearance=0.0,
        )
        hospital.admit(patient, hospital_id=_H1, department=Department.WARD)
        for _ in range(10):
            hospital.step()
        assert patient.state == HealthState.CARRIER
        req = patient.make_micro_request(run_id="run_compat", day=10, dt_days=1, seed=42)
        assert req is not None
        assert req["patient_id"] == "compat_001"
        assert req["episode_id"] == "ep_compat"

    def test_micro_request_valid_after_transfer(self, hospital_network: MacroSimulator):
        """A transferred carrier still produces a valid micro request."""
        patient = Patient(
            patient_id="xfer_compat",
            state=HealthState.CARRIER,
            episode_id="ep_xfer",
            p_clearance=0.0,
        )
        hospital_network.admit(patient, hospital_id=_H1, department=Department.WARD)
        hospital_network.step()
        hospital_network.transfer(patient, to_hospital_id=_H2)
        hospital_network.step()
        req = patient.make_micro_request(run_id="run_xfer", day=2, dt_days=1, seed=42)
        assert req is not None
        assert req["patient_id"] == "xfer_compat"

    def test_micro_request_reflects_icu_after_admission(self, hospital: MacroSimulator):
        """Micro request setting field must reflect the ICU department after admission."""
        patient = Patient(
            patient_id="adm_compat",
            state=HealthState.CARRIER,
            episode_id="ep_adm",
        )
        hospital.admit(patient, hospital_id=_H1, department=Department.ICU)
        hospital.step()
        req = patient.make_micro_request(run_id="run_adm", day=1, dt_days=1, seed=42)
        assert req is not None
        assert req["setting"] == "icu"

    def test_context_update_preserves_micro_fields(self, hospital: MacroSimulator):
        """Repeated update_context (via step) must not overwrite micro-derived fields."""
        patient = Patient(
            patient_id="no_corrupt",
            state=HealthState.CARRIER,
            episode_id="ep_nc",
            resistant_fraction=0.7,
            dominant_genotype="R2",
            relative_transmissibility=1.8,
            p_clearance=0.0,
            severity_modifier=1.5,
            lethality_modifier=1.3,
        )
        hospital.admit(patient, hospital_id=_H1, department=Department.WARD)
        for _ in range(10):
            hospital.step()
        assert patient.resistant_fraction == 0.7
        assert patient.dominant_genotype == "R2"
        assert patient.relative_transmissibility == 1.8
        assert patient.p_clearance == 0.0
        assert patient.severity_modifier == 1.5
        assert patient.lethality_modifier == 1.3

    def test_public_helper_methods(self, carrier_patient: Patient, susceptible_patient: Patient):
        """Patient helper methods used by macro must be callable and well-typed."""
        rng = random.Random(42)
        # transmission_multiplier_for_macro
        tm = carrier_patient.transmission_multiplier_for_macro()
        assert isinstance(tm, float) and tm >= 0.0
        # susceptibility_multiplier_for_macro
        sm = susceptible_patient.susceptibility_multiplier_for_macro()
        assert isinstance(sm, float) and sm >= 0.0
        # should_clear_today
        assert isinstance(carrier_patient.should_clear_today(rng), bool)

    def test_clear_carriage_resets_state(self, carrier_patient: Patient):
        """Patient.clear_carriage() must reset carrier state to susceptible defaults."""
        carrier_patient.clear_carriage()
        assert carrier_patient.state == HealthState.SUSCEPTIBLE
        assert carrier_patient.episode_id is None
        assert carrier_patient.resistant_fraction == 0.0

    def test_discharge_preserves_micro_request_data(self, hospital: MacroSimulator):
        """Discharge + re-admit must not corrupt data required by micro requests."""
        patient = Patient(
            patient_id="dc_compat",
            state=HealthState.CARRIER,
            episode_id="ep_dc",
            resistant_fraction=0.4,
            dominant_genotype="R1",
            p_clearance=0.0,
        )
        hospital.admit(patient, hospital_id=_H1, department=Department.WARD)
        hospital.step()
        hospital.discharge(patient)
        hospital.admit(patient, hospital_id=_H1, department=Department.ICU)
        hospital.step()
        req = patient.make_micro_request(run_id="run_dc", day=1, dt_days=1, seed=42)
        assert req is not None
        assert req["initial_state"]["resistant_fraction"] == 0.4
        assert req["initial_state"]["dominant_genotype"] == "R1"

    @pytest.mark.parametrize("state", [HealthState.SUSCEPTIBLE, HealthState.CARRIER])
    def test_only_s_and_c_states_in_macro(self, state: HealthState, hospital: MacroSimulator):
        """After 10 days every patient must still be S or C — no infection state."""
        episode = "ep_sc" if state == HealthState.CARRIER else None
        patient = Patient(
            patient_id=f"state_{state.value}",
            state=state,
            episode_id=episode,
        )
        hospital.admit(patient, hospital_id=_H1, department=Department.WARD)
        for _ in range(10):
            hospital.step()
        assert patient.state in _VALID_STATES
