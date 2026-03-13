"""
MacroSimulator — hospital-level carrier transmission and clearance.

Responsibilities:
- admit / discharge / transfer patients across hospitals
- daily PatientDailyContext updates (hygiene, isolation, ABX regimen)
- stochastic S->C transmission based on carrier pressure and patient modifiers
- stochastic C->S clearance via Patient.should_clear_today / clear_carriage
- seeded randomness for full reproducibility
"""

from __future__ import annotations

import random
from typing import Dict, List

from exchange.patient import (
    AntibioticRegimen,
    Department,
    HealthState,
    Patient,
    PatientDailyContext,
)
from macro_simulation.simulation import SimulationConfig

# Antibiotic classes the macro policy may prescribe
_ABX_CLASSES = [
    "beta_lactam",
    "fluoroquinolone",
    "glycopeptide",
    "macrolide",
    "aminoglycoside",
]
_DOSE_LEVELS = ["low", "std", "high"]


class MacroSimulator:
    """Simulates hospital-level transmission across one or more hospitals.

    Parameters
    ----------
    config : SimulationConfig
        Simulation hyper-parameters (hygiene, ABX policy, …).
    n_hospitals : int
        Number of hospitals in the network.
    seed : int | None
        RNG seed for reproducibility.
    """

    def __init__(
        self,
        config: SimulationConfig,
        n_hospitals: int = 1,
        seed: int | None = None,
    ) -> None:
        self._config = config
        self._rng = random.Random(seed)

        # Build hospital registry
        self._hospital_ids = [f"hospital_{i:03d}" for i in range(1, n_hospitals + 1)]
        # hospital_id -> list of Patient objects currently admitted
        self._patients: Dict[str, List[Patient]] = {hid: [] for hid in self._hospital_ids}
        # patient_id -> hospital_id (fast lookup)
        self._patient_hospital: Dict[str, str] = {}
        # patient_id -> department
        self._patient_department: Dict[str, Department] = {}

        self._day = 0
        self._episode_counter = 0

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def admit(
        self,
        patient: Patient,
        hospital_id: str,
        department: Department,
    ) -> None:
        """Admit patient to hospital_id in the given department."""
        if hospital_id not in self._patients:
            raise ValueError(f"Unknown hospital: {hospital_id}")

        patient.hospital_id = hospital_id
        patient.department = department

        self._patients[hospital_id].append(patient)
        self._patient_hospital[patient.patient_id] = hospital_id
        self._patient_department[patient.patient_id] = department

    def discharge(self, patient: Patient) -> None:
        """Remove patient from its current hospital."""
        hid = self._patient_hospital.get(patient.patient_id)
        if hid is None:
            raise ValueError(f"Patient {patient.patient_id} is not currently admitted.")

        self._patients[hid].remove(patient)
        del self._patient_hospital[patient.patient_id]
        del self._patient_department[patient.patient_id]
        patient.hospital_id = None

    def transfer(self, patient: Patient, to_hospital_id: str) -> None:
        """Transfer patient to to_hospital_id."""
        src = self._patient_hospital.get(patient.patient_id)
        if src is None:
            raise ValueError(f"Patient {patient.patient_id} is not currently admitted.")
        if to_hospital_id == src:
            raise ValueError("Cannot transfer to the same hospital.")
        if to_hospital_id not in self._patients:
            raise ValueError(f"Unknown hospital: {to_hospital_id}")

        self._patients[src].remove(patient)
        self._patients[to_hospital_id].append(patient)
        self._patient_hospital[patient.patient_id] = to_hospital_id
        patient.hospital_id = to_hospital_id

    def get_occupancy(self, hospital_id: str) -> int:
        """Return the number of patients currently in hospital_id."""
        return len(self._patients.get(hospital_id, []))

    def get_patients(self, hospital_id: str) -> List[Patient]:
        """Return a snapshot of patients currently in hospital_id."""
        return list(self._patients.get(hospital_id, []))

    # ------------------------------------------------------------------
    # Daily step
    # ------------------------------------------------------------------

    def step(self) -> None:
        """Advance the simulation by one day.

        Order of operations per hospital:
        1. Clearance: each carrier may revert to susceptible (C -> S).
        2. Update PatientDailyContext for every patient.
        3. Transmission: susceptible patients may become carriers (S -> C).
        """
        self._day += 1

        for hid in self._hospital_ids:
            patients = self._patients[hid]
            if not patients:
                continue

            # 1. Clearance (C -> S)
            for p in patients:
                if p.state == HealthState.CARRIER:
                    if p.should_clear_today(self._rng):
                        p.clear_carriage()

            # 2. Context update
            for p in patients:
                ctx = self._build_context(p, hid)
                p.update_context(ctx)

            # 3. Transmission (S -> C)
            self._do_transmission(patients)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_context(self, patient: Patient, hospital_id: str) -> PatientDailyContext:
        """Construct the daily context snapshot for patient.

        Always consumes exactly 4 random draws from ``self._rng`` regardless
        of branching, so that downstream draws (clearance, transmission) land
        on deterministic positions in the RNG stream.
        """
        department = self._patient_department[patient.patient_id]
        cfg = self._config

        # Draw all random values up-front (fixed count = 4)
        iso_draw = self._rng.random()
        abx_draw = self._rng.random()
        class_draw = self._rng.random()
        dose_draw = self._rng.random()

        # --- Isolation decision ---
        is_isolated = patient.is_isolated
        if not is_isolated and patient.state == HealthState.CARRIER:
            if iso_draw < cfg.carrier_isolation_probability:
                is_isolated = True

        # --- Antibiotic regimen ---
        abx_prob = (
            cfg.icu_abx_probability if department == Department.ICU else cfg.ward_abx_probability
        )
        if abx_draw < abx_prob:
            class_idx = int(class_draw * len(_ABX_CLASSES)) % len(_ABX_CLASSES)
            dose_idx = int(dose_draw * len(_DOSE_LEVELS)) % len(_DOSE_LEVELS)
            regimen = AntibioticRegimen(
                on=True,
                abx_class=_ABX_CLASSES[class_idx],
                dose_level=_DOSE_LEVELS[dose_idx],
            )
        else:
            regimen = AntibioticRegimen(on=False, abx_class="none", dose_level="std")

        return PatientDailyContext(
            hospital_id=hospital_id,
            department=department,
            hygiene_level=cfg.base_hygiene,
            isolation_effectiveness=cfg.base_isolation_effectiveness,
            diagnostic_speed=cfg.base_diagnostic_speed,
            is_isolated=is_isolated,
            regimen=regimen,
        )

    def _do_transmission(self, patients: List[Patient]) -> None:
        """Stochastic S -> C transmission within one hospital ward."""
        carriers = [p for p in patients if p.state == HealthState.CARRIER]
        susceptible = [p for p in patients if p.state == HealthState.SUSCEPTIBLE]

        if not carriers or not susceptible:
            return

        cfg = self._config
        hygiene_factor = 1.0 - cfg.base_hygiene  # higher hygiene => lower factor
        iso_reduction = cfg.base_isolation_effectiveness

        # Sum up carrier-side force of infection
        carrier_force = 0.0
        for car in carriers:
            contribution = cfg.base_transmission_rate * car.transmission_multiplier_for_macro()
            if car.is_isolated:
                contribution *= 1.0 - iso_reduction
            carrier_force += contribution

        carrier_force *= hygiene_factor

        # For each susceptible patient, draw colonisation
        for sus in susceptible:
            p_colonize = carrier_force * sus.susceptibility_multiplier_for_macro()
            if sus.is_isolated:
                p_colonize *= 1.0 - iso_reduction
            p_colonize = min(1.0, p_colonize)

            if self._rng.random() < p_colonize:
                self._colonize(sus)

    def _colonize(self, patient: Patient) -> None:
        """Transition a susceptible patient to carrier (S -> C)."""
        self._episode_counter += 1
        patient.state = HealthState.CARRIER
        patient.episode_id = f"episode_new_{self._episode_counter}"
