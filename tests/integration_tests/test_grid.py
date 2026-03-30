from __future__ import annotations


from exchange.patient import Department, HealthState, Patient
from macro_simulation.agents import PatientAgent, _MesaModel
from macro_simulation.grid import HospitalDepartmentGrid, HospitalNetworkGrid
from macro_simulation.simulation import SimulationConfig
from macro_simulation.simulator import MacroSimulator

_H1 = "hospital_001"
_H2 = "hospital_002"


class _FixedChoiceRng:
    def choice(self, seq):
        return seq[0]


def test_assign_ward_cell_in_ward_zone() -> None:
    model = _MesaModel()
    grid = HospitalDepartmentGrid(
        cols=3,
        rows=3,
        icu_rows=1,
        mesa_model=model,
        rng=_FixedChoiceRng(),
        alpha=0.5,
    )
    agent = PatientAgent(Patient(patient_id="p_ward"), model)
    pos = grid.assign_to_department(agent, Department.WARD)

    assert grid.dept_type_at(pos) == Department.WARD


def test_assign_icu_cell_in_icu_zone() -> None:
    model = _MesaModel()
    grid = HospitalDepartmentGrid(
        cols=3,
        rows=3,
        icu_rows=1,
        mesa_model=model,
        rng=_FixedChoiceRng(),
        alpha=0.5,
    )
    agent = PatientAgent(Patient(patient_id="p_icu"), model)
    pos = grid.assign_to_department(agent, Department.ICU)

    assert grid.dept_type_at(pos) == Department.ICU


def test_multipatient_same_department_cell() -> None:
    model = _MesaModel()
    grid = HospitalDepartmentGrid(
        cols=2,
        rows=2,
        icu_rows=1,
        mesa_model=model,
        rng=_FixedChoiceRng(),
        alpha=0.5,
    )
    agent_a = PatientAgent(Patient(patient_id="p_a"), model)
    agent_b = PatientAgent(Patient(patient_id="p_b"), model)

    pos = grid.assign_to_department(agent_a, Department.WARD)
    grid.assign_to_department(agent_b, Department.WARD)

    same_cell = [a for a in grid.get_all_agents() if a.pos == pos]
    assert len(same_cell) == 2


def test_release_removes_patient() -> None:
    model = _MesaModel()
    grid = HospitalDepartmentGrid(
        cols=2,
        rows=2,
        icu_rows=1,
        mesa_model=model,
        rng=_FixedChoiceRng(),
        alpha=0.5,
    )
    agent = PatientAgent(Patient(patient_id="p_release"), model)
    grid.assign_to_department(agent, Department.WARD)

    grid.release(agent)

    assert grid.get_all_patients() == []


def test_proximity_weight_same_cell_is_one() -> None:
    model = _MesaModel()
    grid = HospitalDepartmentGrid(
        cols=2,
        rows=2,
        icu_rows=1,
        mesa_model=model,
        rng=_FixedChoiceRng(),
        alpha=0.5,
    )
    pos = (0, 1)

    assert grid.proximity_weight(pos, pos) == 1.0


def test_proximity_weight_decays_with_chebyshev() -> None:
    model = _MesaModel()
    grid = HospitalDepartmentGrid(
        cols=3,
        rows=2,
        icu_rows=1,
        mesa_model=model,
        rng=_FixedChoiceRng(),
        alpha=0.5,
    )
    w_near = grid.proximity_weight((0, 1), (1, 1))
    w_far = grid.proximity_weight((0, 1), (2, 1))

    assert w_near > w_far


def test_dept_type_at_icu_boundary() -> None:
    model = _MesaModel()
    grid = HospitalDepartmentGrid(
        cols=2,
        rows=2,
        icu_rows=1,
        mesa_model=model,
        rng=_FixedChoiceRng(),
        alpha=0.5,
    )

    assert grid.dept_type_at((0, 0)) == Department.ICU
    assert grid.dept_type_at((0, 1)) == Department.WARD


def test_alpha_zero_approximates_well_mixed() -> None:
    model = _MesaModel()
    grid = HospitalDepartmentGrid(
        cols=3,
        rows=2,
        icu_rows=1,
        mesa_model=model,
        rng=_FixedChoiceRng(),
        alpha=0.0,
    )

    assert grid.proximity_weight((0, 1), (2, 1)) == 1.0


def test_six_hospitals_3x2_positions() -> None:
    ids = [f"hospital_{i:03d}" for i in range(1, 7)]
    grid = HospitalNetworkGrid(ids, cols=3)

    assert grid.get_position("hospital_001") == (0, 0)
    assert grid.get_position("hospital_002") == (1, 0)
    assert grid.get_position("hospital_003") == (2, 0)
    assert grid.get_position("hospital_004") == (0, 1)
    assert grid.get_position("hospital_005") == (1, 1)
    assert grid.get_position("hospital_006") == (2, 1)


def test_adjacent_hospitals_are_neighbors() -> None:
    ids = [f"hospital_{i:03d}" for i in range(1, 7)]
    grid = HospitalNetworkGrid(ids, cols=3)

    neighbors = grid.get_neighbors("hospital_001")
    assert "hospital_002" in neighbors
    assert "hospital_004" in neighbors


def test_admit_sets_agent_pos_and_dept_type() -> None:
    sim = MacroSimulator(config=SimulationConfig(), n_hospitals=1, seed=10)
    patient = Patient(patient_id="p_admit")
    sim.admit(patient, hospital_id=_H1, department=Department.WARD)

    agent = sim._patient_agents[patient.patient_id]
    grid = sim._dept_grids[_H1]

    assert agent.pos is not None
    assert grid.dept_type_at(agent.pos) == Department.WARD


def test_discharge_removes_agent_from_grid() -> None:
    sim = MacroSimulator(config=SimulationConfig(), n_hospitals=1, seed=11)
    patient = Patient(patient_id="p_discharge")
    sim.admit(patient, hospital_id=_H1, department=Department.WARD)

    sim.discharge(patient)

    assert patient.patient_id not in sim._patient_agents
    assert patient not in sim.get_patients(_H1)


def test_transfer_reassigns_to_new_grid() -> None:
    sim = MacroSimulator(config=SimulationConfig(), n_hospitals=2, seed=12)
    patient = Patient(patient_id="p_transfer")
    sim.admit(patient, hospital_id=_H1, department=Department.WARD)

    sim.transfer(patient, to_hospital_id=_H2)

    assert patient.hospital_id == _H2
    assert patient in sim.get_patients(_H2)
    assert patient not in sim.get_patients(_H1)


def test_get_patients_uses_grid() -> None:
    sim = MacroSimulator(config=SimulationConfig(), n_hospitals=1, seed=13)
    patient = Patient(patient_id="p_grid")
    sim.admit(patient, hospital_id=_H1, department=Department.WARD)

    assert sim.get_patients(_H1) == [patient]


def test_proximity_transmission_higher_same_dept() -> None:
    cfg = SimulationConfig(
        base_hygiene=0.0,
        base_transmission_rate=0.4,
        daily_contact_attempts=1.0,
        proximity_decay_alpha=2.0,
        dept_grid_cols=3,
        dept_grid_rows=2,
        dept_grid_icu_rows=1,
    )

    def run_trial(seed: int, same_cell: bool) -> bool:
        sim = MacroSimulator(config=cfg, n_hospitals=1, seed=seed)
        carrier = Patient(
            patient_id=f"car_{seed}",
            state=HealthState.CARRIER,
            episode_id=f"ep_{seed}",
            p_clearance=0.0,
        )
        susceptible = Patient(patient_id=f"sus_{seed}", state=HealthState.SUSCEPTIBLE)
        sim.admit(carrier, hospital_id=_H1, department=Department.WARD)
        sim.admit(susceptible, hospital_id=_H1, department=Department.WARD)

        grid = sim._dept_grids[_H1]
        car_agent = sim._patient_agents[carrier.patient_id]
        sus_agent = sim._patient_agents[susceptible.patient_id]
        if same_cell:
            target = (0, 1)
            grid._grid.move_agent(car_agent, target)
            grid._grid.move_agent(sus_agent, target)
        else:
            grid._grid.move_agent(car_agent, (0, 1))
            grid._grid.move_agent(sus_agent, (2, 1))

        sim._do_transmission(_H1)
        return susceptible.state == HealthState.CARRIER

    same_count = sum(run_trial(seed, True) for seed in range(200))
    far_count = sum(run_trial(seed + 1000, False) for seed in range(200))

    assert same_count >= far_count + 10
