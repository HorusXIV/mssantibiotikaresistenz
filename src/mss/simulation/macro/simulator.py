"""
MacroSimulator — hospital-level carrier transmission and clearance.

Responsibilities:
- admit / discharge / transfer patients across hospitals
- daily PatientDailyContext updates (hygiene, isolation, ABX regimen)
- optional daily micro exchange for all active carrier episodes
- stochastic S->C transmission based on carrier pressure and patient modifiers
- stochastic C->S clearance via Patient.should_clear_today / clear_carriage
- seeded randomness for full reproducibility
"""

from __future__ import annotations

import math
import random
from typing import Any, Callable, Dict, List, Optional

from mss.domain import (
    AntibioticRegimen,
    Department,
    HealthState,
    Patient,
    PatientDailyContext,
)
from mss.simulation.macro.agents import PatientAgent, _MesaModel
from mss.simulation.macro.config import SimulationConfig
from mss.simulation.macro.grid import HospitalDepartmentGrid, HospitalNetworkGrid

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
        # patient_id -> hospital_id (fast lookup)
        self._patient_hospital: Dict[str, str] = {}

        base_seed = seed if seed is not None else None
        self._mesa_model = _MesaModel()
        self._dept_grids: Dict[str, HospitalDepartmentGrid] = {
            hid: HospitalDepartmentGrid(
                cols=config.dept_grid_cols,
                rows=config.dept_grid_rows,
                icu_rows=config.dept_grid_icu_rows,
                mesa_model=self._mesa_model,
                rng=random.Random(None if base_seed is None else base_seed + 1000 + i),
                alpha=config.proximity_decay_alpha,
            )
            for i, hid in enumerate(self._hospital_ids)
        }
        self._network_grid = HospitalNetworkGrid(
            self._hospital_ids,
            cols=config.network_grid_cols,
            mesa_model=self._mesa_model,
        )
        self._patient_agents: Dict[str, PatientAgent] = {}

        self._day = 0
        self._episode_counter = 0
        self._daily_transfers: List[Dict[str, str]] = []
        # In-hospital transmission counters for the spatial calibration:
        # how many S->C events had their source in the SAME grid cell (roommate)
        # vs. total. Their ratio is set by proximity_decay_alpha.
        self._transmission_count = 0
        self._same_cell_transmission_count = 0

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
        if hospital_id not in self._dept_grids:
            raise ValueError(f"Unknown hospital: {hospital_id}")

        patient.hospital_id = hospital_id
        patient.department = department

        self._patient_hospital[patient.patient_id] = hospital_id
        agent = PatientAgent(patient, self._mesa_model)
        self._patient_agents[patient.patient_id] = agent
        self._dept_grids[hospital_id].assign_to_department(agent, department)

        patient.admission_day = self._day
        patient.planned_discharge_day = self._day + self._draw_los(department)

    def discharge(self, patient: Patient, micro_simulator: Any = None) -> None:
        """Remove patient from its current hospital.

        If ``micro_simulator`` is provided and the patient is a carrier, the
        associated micro episode is cleaned up before removal.
        """
        hid = self._patient_hospital.get(patient.patient_id)
        if hid is None:
            raise ValueError(f"Patient {patient.patient_id} is not currently admitted.")

        if micro_simulator is not None and patient.state == HealthState.CARRIER:
            if patient.episode_id is not None:
                clear_fn = getattr(micro_simulator, "clear_episode", None)
                if clear_fn is not None:
                    clear_fn(patient.episode_id)

        agent = self._patient_agents.pop(patient.patient_id, None)
        if agent is not None:
            self._dept_grids[hid].release(agent)
        del self._patient_hospital[patient.patient_id]
        patient.hospital_id = None

    def transfer(self, patient: Patient, to_hospital_id: str) -> None:
        """Transfer patient to to_hospital_id."""
        src = self._patient_hospital.get(patient.patient_id)
        if src is None:
            raise ValueError(f"Patient {patient.patient_id} is not currently admitted.")
        if to_hospital_id == src:
            raise ValueError("Cannot transfer to the same hospital.")
        if to_hospital_id not in self._dept_grids:
            raise ValueError(f"Unknown hospital: {to_hospital_id}")

        agent = self._patient_agents.get(patient.patient_id)
        if agent is None:
            raise ValueError(f"Patient {patient.patient_id} is missing a grid agent.")

        src_grid = self._dept_grids[src]
        dst_grid = self._dept_grids[to_hospital_id]
        current_dept = src_grid.dept_type_at(agent.pos)
        src_grid.release(agent)
        dst_grid.assign_to_department(agent, current_dept)

        self._patient_hospital[patient.patient_id] = to_hospital_id
        patient.hospital_id = to_hospital_id
        patient.department = current_dept

    def get_occupancy(self, hospital_id: str) -> int:
        """Return the number of patients currently in hospital_id."""
        grid = self._dept_grids.get(hospital_id)
        if grid is None:
            return 0
        return len(grid.get_all_patients())

    def get_patients(self, hospital_id: str) -> List[Patient]:
        """Return a snapshot of patients currently in hospital_id."""
        grid = self._dept_grids.get(hospital_id)
        if grid is None:
            return []
        return list(grid.get_all_patients())

    def get_cell_stats(self, hospital_id: str) -> List[Dict[str, object]]:
        """Per-cell occupancy (x, y, department, total, carriers, prevalence)."""
        grid = self._dept_grids.get(hospital_id)
        if grid is None:
            return []
        return grid.cell_stats()

    def get_colonization_count(self) -> int:
        """Cumulative count of in-hospital S->C transmissions (true nosocomial
        acquisitions). Incremented only in _colonize(); community/seed carriers
        admitted directly as carriers do not count. The per-day delta is the
        clean acquisition incidence, free of community-import and discharge churn.
        """
        return self._episode_counter

    def get_transmission_count(self) -> int:
        """Cumulative number of in-hospital transmission events (S->C)."""
        return self._transmission_count

    def get_same_cell_transmission_count(self) -> int:
        """Cumulative in-hospital transmissions whose source was in the same grid
        cell (roommate). same_cell / total is set by proximity_decay_alpha and is
        the spatial observable to calibrate it against (roommate-attributable
        transmission ~20-30% in the literature)."""
        return self._same_cell_transmission_count

    # ------------------------------------------------------------------
    # Daily step
    # ------------------------------------------------------------------

    def step(
        self,
        micro_simulator: Any = None,
        run_id: str = "run_001",
        patient_factory: Optional[Callable[[str, Department], Patient]] = None,
    ) -> None:
        """Advance the simulation by one day.

        Order of operations per day:
        1. Per hospital clearance: each carrier may revert to susceptible (C -> S).
        2. Per hospital LOS-based discharge (logistic probability after planned date).
           Clearance runs first so a carrier who clears today and is already past
           their planned discharge date can leave the same day.
        3. Poisson admissions to least-occupied hospital (capped).
        4. Per hospital PatientDailyContext update for every patient.
           Newly detected carriers (is_isolated set True) immediately have their
           planned_discharge_day extended by carrier_extension_days.
        5. One optional global micro batch for all active carriers.
        6. Per hospital transmission: susceptible patients may become carriers (S -> C).

        Parameters
        ----------
        patient_factory
            Optional callable ``(hospital_id, department) -> Patient`` used by
            the admission phase. If None or ``daily_admission_rate == 0``,
            the admission phase is skipped.
        """
        self._day += 1

        all_micro_requests: List[Dict[str, Any]] = []
        request_patients: List[Patient] = []
        active_hospital_ids: List[str] = []

        for hid in self._hospital_ids:
            patients = self._dept_grids[hid].get_all_patients()
            if not patients:
                continue
            active_hospital_ids.append(hid)

            # 1a. Clearance (C -> S): must run before discharge so a carrier
            #     who clears today is seen as susceptible in the discharge phase.
            for p in patients:
                if p.state == HealthState.CARRIER:
                    if p.should_clear_today(self._rng):
                        if micro_simulator is not None and p.episode_id is not None:
                            clear_episode = getattr(micro_simulator, "clear_episode", None)
                            if clear_episode is not None:
                                clear_episode(p.episode_id)
                        p.clear_carriage()

            # 1b. Context update + isolation detection
            for p in patients:
                was_isolated = p.is_isolated
                ctx = self._build_context(p, hid)
                p.update_context(ctx)
                # On first detection, push discharge to detection day +
                # carrier_extension_days (scaled by severity_modifier). The clock
                # runs from detection and is not rolled forward again.
                if p.state == HealthState.CARRIER and p.is_isolated and not was_isolated:
                    extension = max(
                        1, int(self._config.carrier_extension_days * p.severity_modifier)
                    )
                    extended = self._day + extension
                    if p.planned_discharge_day is None or extended > p.planned_discharge_day:
                        p.planned_discharge_day = extended

            # 1c. Collect micro requests (carrier episodes)
            if micro_simulator is not None:
                requests, batch_patients = self._collect_micro_requests(
                    patients=patients, run_id=run_id
                )
                all_micro_requests.extend(requests)
                request_patients.extend(batch_patients)

        # 2. Discharge (after clearance, so post-clearance state is visible here)
        for hid in self._hospital_ids:
            self._discharge_phase(
                list(self._dept_grids[hid].get_all_patients()), hid, micro_simulator
            )

        # 3. Admissions
        self._admission_phase(patient_factory)

        # 4. Micro batch
        if micro_simulator is not None and all_micro_requests:
            # Requests were collected before the discharge phase. Drop any whose
            # patient left today; discharge() already cleared their episode, so
            # re-running the micro batch would re-store it and leak episode state.
            kept = [
                (req, pat)
                for req, pat in zip(all_micro_requests, request_patients)
                if pat.hospital_id is not None
            ]
            if kept:
                reqs, pats = map(list, zip(*kept))
                self._apply_micro_phase(
                    micro_simulator=micro_simulator,
                    requests=reqs,
                    request_patients=pats,
                )

        # 5. Transmission (S -> C)
        for hid in active_hospital_ids:
            self._do_transmission(hid)

        # 6. Inter-hospital transfers
        self._transfer_phase()

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _draw_los(self, department: Department) -> int:
        """Draw a length-of-stay from a log-normal distribution."""
        cfg = self._config
        mean = cfg.los_mean_icu if department == Department.ICU else cfg.los_mean_ward
        sigma = cfg.los_sigma
        mu = math.log(mean) - 0.5 * sigma**2
        return max(1, round(self._rng.lognormvariate(mu, sigma)))

    def _discharge_phase(
        self,
        patients: List[Patient],
        hospital_id: str,
        micro_simulator: Any,
    ) -> None:
        """Discharge eligible patients via the logistic length-of-stay curve.

        All patients face a daily mortality probability (base_mortality_rate);
        for carriers this is scaled by severity_modifier — more virulent strains
        carry higher mortality risk. Carriers and susceptibles share the same
        logistic discharge once they pass their planned discharge day. Detected
        carriers already had that day pushed back by carrier_extension_days (in
        the context phase), so they stay longer but are still discharged — they
        leave colonized, mirroring real MRSA carriers sent home under contact
        precautions.
        """
        cfg = self._config
        admissions_active = cfg.daily_admission_rate > 0.0
        for p in patients:
            if p.planned_discharge_day is None:
                continue

            # Mortality is evaluated daily, independent of admissions. Skipping the
            # draw at rate 0 keeps the RNG stream stable for closed-cohort runs.
            mortality_rate = cfg.base_mortality_rate
            if p.state == HealthState.CARRIER:
                mortality_rate *= p.severity_modifier
            if mortality_rate > 0.0 and self._rng.random() < mortality_rate:
                self.discharge(p, micro_simulator)
                continue

            # Logistic length-of-stay discharge runs only while admissions are
            # active (a closed cohort keeps everyone for clean observation).
            if not admissions_active:
                continue

            if self._day >= p.planned_discharge_day:
                days_over = self._day - p.planned_discharge_day
                p_d = 1.0 / (
                    1.0
                    + math.exp(
                        -cfg.discharge_logistic_k * (days_over - cfg.discharge_logistic_t_half)
                    )
                )
                if self._rng.random() < p_d:
                    self.discharge(p, micro_simulator)

    @staticmethod
    def _poisson_draw(rng: random.Random, lam: float) -> int:
        """Knuth's algorithm for Poisson sampling without numpy.

        lam is capped at 700 so exp(-lam) cannot underflow to 0.0 (which would
        loop forever) for any realistic admission rate.
        """
        if lam <= 0:
            return 0
        big_l = math.exp(-min(lam, 700.0))
        k, p = 0, 1.0
        while p > big_l:
            k += 1
            p *= rng.random()
        return k - 1

    def _admission_phase(
        self,
        patient_factory: Optional[Callable[[str, Department], Patient]],
    ) -> None:
        """Admit Poisson-drawn new patients to the least-occupied hospital under cap."""
        cfg = self._config
        if patient_factory is None or cfg.daily_admission_rate <= 0.0:
            return
        cap = cfg.max_occupancy_per_hospital
        n_new = self._poisson_draw(self._rng, cfg.daily_admission_rate)
        for _ in range(n_new):
            candidates = [h for h in self._hospital_ids if self.get_occupancy(h) < cap]
            if not candidates:
                break
            # Weight by remaining free capacity — less full = more likely, but stochastic
            weights = [cap - self.get_occupancy(h) for h in candidates]
            hid = self._rng.choices(candidates, weights=weights, k=1)[0]
            dept = self._draw_admission_department(self._dept_grids[hid])
            new_patient = patient_factory(hid, dept)
            self.admit(new_patient, hid, dept)

    def _draw_admission_department(self, grid: HospitalDepartmentGrid) -> Department:
        ward_positions = grid.dept_positions_for(Department.WARD)
        icu_positions = grid.dept_positions_for(Department.ICU)
        ward_weight = len(ward_positions)
        icu_weight = len(icu_positions)

        if ward_weight <= 0 and icu_weight <= 0:
            return Department.WARD
        if icu_weight <= 0:
            return Department.WARD
        if ward_weight <= 0:
            return Department.ICU

        if self._rng.random() < (icu_weight / (icu_weight + ward_weight)):
            return Department.ICU
        return Department.WARD

    def _build_context(self, patient: Patient, hospital_id: str) -> PatientDailyContext:
        """Construct the daily context snapshot for patient.

        Always consumes exactly 4 random draws from ``self._rng`` regardless
        of branching, so that downstream draws (clearance, transmission) land
        on deterministic positions in the RNG stream.
        """
        agent = self._patient_agents.get(patient.patient_id)
        if agent is None:
            raise ValueError(f"Patient {patient.patient_id} is missing a grid agent.")
        department = self._dept_grids[hospital_id].dept_type_at(agent.pos)
        cfg = self._config

        # Draw all random values up-front (fixed count = 4)
        iso_draw = self._rng.random()
        abx_draw = self._rng.random()
        class_draw = self._rng.random()
        dose_draw = self._rng.random()

        # --- Isolation decision ---
        # diagnostic_speed scales the effective daily detection probability:
        # 1.0 (default) = unchanged, 0.5 = half speed, 2.0 = double (capped at 1.0)
        effective_detection_prob = min(
            1.0, cfg.carrier_isolation_probability * cfg.base_diagnostic_speed
        )
        is_isolated = patient.is_isolated
        if not is_isolated and patient.state == HealthState.CARRIER:
            if iso_draw < effective_detection_prob:
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

    def _do_transmission(self, hospital_id: str) -> None:
        """Stochastic S -> C transmission within one hospital (proximity-weighted)."""
        grid = self._dept_grids[hospital_id]
        agents = grid.get_all_agents()
        carriers = [a for a in agents if a.patient.state == HealthState.CARRIER]
        susceptible = [a for a in agents if a.patient.state == HealthState.SUSCEPTIBLE]

        if not carriers or not susceptible:
            return

        # Source position per carrier (for the same-cell / roommate transmission metric).
        carrier_pos = {a.patient.patient_id: a.pos for a in carriers}

        cfg = self._config
        n_total = max(1, len(agents))

        for sus_agent in susceptible:
            sus_pos = sus_agent.pos
            sus = sus_agent.patient

            # Read hygiene from the susceptible patient's daily context (set by _build_context).
            # Falls back to the global config value if context is not yet initialised.
            sus_ctx = sus._ctx
            hygiene = sus_ctx.hygiene_level if sus_ctx is not None else cfg.base_hygiene
            hygiene_factor = 1.0 - hygiene

            weighted_force = 0.0
            weighted_carriers: List[tuple[Patient, float]] = []
            for car_agent in carriers:
                w = grid.proximity_weight(sus_pos, car_agent.pos)
                contrib = (
                    cfg.base_transmission_rate
                    * car_agent.patient.transmission_multiplier_for_macro()
                    * w
                )
                if car_agent.patient.is_isolated:
                    # Use the carrier's own isolation effectiveness from their daily context.
                    car_ctx = car_agent.patient._ctx
                    iso_eff = (
                        car_ctx.isolation_effectiveness
                        if car_ctx is not None
                        else cfg.base_isolation_effectiveness
                    )
                    contrib *= 1.0 - iso_eff
                weighted_force += contrib
                weighted_carriers.append((car_agent.patient, contrib))

            # Daily infection hazard: contacts x (1 - hygiene) x normalized carrier
            # pressure x host susceptibility; p_colonize = 1 - exp(-hazard).
            hazard = (
                cfg.daily_contact_attempts
                * hygiene_factor
                * (weighted_force / n_total)
                * sus.susceptibility_multiplier_for_macro()
            )
            if sus.is_isolated:
                sus_iso_eff = (
                    sus_ctx.isolation_effectiveness
                    if sus_ctx is not None
                    else cfg.base_isolation_effectiveness
                )
                hazard *= 1.0 - sus_iso_eff

            p_colonize = 1.0 - math.exp(-max(0.0, hazard))
            if self._rng.random() < p_colonize:
                source = self._pick_transmission_source(weighted_carriers, weighted_force)
                self._transmission_count += 1
                if source is not None:
                    src_pos = carrier_pos.get(source.patient_id)
                    if (
                        src_pos is not None
                        and grid.chebyshev(sus_pos[0], sus_pos[1], src_pos[0], src_pos[1]) == 0
                    ):
                        self._same_cell_transmission_count += 1
                self._colonize(sus, source=source)

    def get_daily_transfers(self) -> List[Dict[str, str]]:
        """Return the list of transfers that occurred during the last step."""
        return list(self._daily_transfers)

    def _transfer_phase(self) -> None:
        """Transfer patients between hospitals, weighted by inverse distance and free capacity."""
        self._daily_transfers = []
        rate = self._config.daily_transfer_rate
        if rate <= 0.0:
            return

        cap = self._config.max_occupancy_per_hospital
        transferred_ids: set[str] = set()

        for hid in self._hospital_ids:
            destinations = [h for h in self._hospital_ids if h != hid]
            if not destinations:
                continue
            for patient in list(self._dept_grids[hid].get_all_patients()):
                if patient.patient_id in transferred_ids:
                    continue
                if self._rng.random() >= rate:
                    continue
                candidates = [d for d in destinations if self.get_occupancy(d) < cap]
                if not candidates:
                    continue
                weights = [
                    (cap - self.get_occupancy(d)) / self._network_grid.distance(hid, d)
                    for d in candidates
                ]
                dest = self._rng.choices(candidates, weights=weights, k=1)[0]
                self.transfer(patient, dest)
                self._daily_transfers.append({"from_hospital": hid, "to_hospital": dest})
                transferred_ids.add(patient.patient_id)

    def _collect_micro_requests(
        self, patients: List[Patient], run_id: str
    ) -> tuple[List[Dict[str, Any]], List[Patient]]:
        """Build one daily micro request per active carrier."""
        requests: List[Dict[str, Any]] = []
        request_patients: List[Patient] = []

        for patient in patients:
            req = patient.make_micro_request(
                run_id=run_id,
                day=self._day,
                dt_days=1,
                seed=self._rng.randint(0, 2**31 - 1),
            )
            if req is not None:
                requests.append(req)
                request_patients.append(patient)

        return requests, request_patients

    def _apply_micro_phase(
        self,
        micro_simulator: Any,
        requests: List[Dict[str, Any]],
        request_patients: List[Patient],
    ) -> None:
        """Run one global micro batch for all active carrier episodes this day."""
        if not requests:
            return

        responses = micro_simulator.process_batch(requests, parallel=True)
        if len(responses) != len(request_patients):
            raise RuntimeError(
                "Micro simulator returned a different number of responses than requests."
            )

        for patient, response in zip(request_patients, responses):
            patient.apply_micro_response(response)

    def _pick_transmission_source(
        self,
        weighted_carriers: List[tuple[Patient, float]],
        total_force: float,
    ) -> Patient | None:
        """Sample which carrier seeded a new colonization event."""
        if total_force <= 0.0:
            return None

        threshold = self._rng.random() * total_force
        running = 0.0
        for carrier, weight in weighted_carriers:
            if weight <= 0.0:
                continue
            running += weight
            if running >= threshold:
                return carrier

        return weighted_carriers[-1][0] if weighted_carriers else None

    def _inherit_transmitted_state(self, source: Patient) -> tuple[float, str]:
        """Copy the source strain to the newly colonized patient."""
        resistance = min(1.0, max(0.0, source.resistant_fraction))
        genotype = source.dominant_genotype
        return resistance, genotype

    def _colonize(self, patient: Patient, source: Patient | None = None) -> None:
        """Transition a susceptible patient to carrier (S -> C)."""
        self._episode_counter += 1
        patient.state = HealthState.CARRIER
        patient.episode_id = f"episode_new_{self._episode_counter}"
        patient.is_isolated = False

        if source is None:
            patient.resistant_fraction = 0.0
            patient.dominant_genotype = "S"
            patient.dominant_strain_name = ""
            patient.relative_transmissibility = 1.0
            patient.severity_modifier = 1.0
            patient.p_clearance = 0.02
            return

        resistance, genotype = self._inherit_transmitted_state(source)
        patient.resistant_fraction = resistance
        patient.dominant_genotype = genotype
        patient.dominant_strain_name = source.dominant_strain_name
        patient.dominant_genome = source.dominant_genome
        patient.relative_transmissibility = source.relative_transmissibility
        patient.severity_modifier = source.severity_modifier
        # Clearance is host-specific and will be recomputed by micro on the next day.
        patient.p_clearance = 0.02
