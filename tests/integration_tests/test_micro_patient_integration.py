"""
Integration tests for Patient <-> MicroSimulator interface compatibility.

These tests verify that:
1. Patient.make_micro_request() produces valid requests for MicroSimulator
2. MicroSimulator responses are compatible with Patient.apply_micro_response()
3. The full data flow works end-to-end without breaking changes
"""

from __future__ import annotations

import numpy as np
import pytest

from mss.domain import (
    AntibioticRegimen,
    Department,
    HealthState,
    Patient,
    PatientDailyContext,
)
from mss.simulation.micro import (
    GeneIndex,
    MicroSimulator,
    SimulationConfig,
    StrainPopulation,
    build_founder_pool,
    create_wild_type_genome,
    get_dominant_strain,
    simulate_day,
)

# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def carrier_patient() -> Patient:
    """Create a patient in CARRIER state with all required fields."""
    patient = Patient(
        patient_id="test_patient_001",
        state=HealthState.CARRIER,
        episode_id="episode_001",
        age_years=55,
        compliance=0.85,
        vulnerability=1.2,
        immune_strength=0.9,
        immune_status="normal",
        sociability=1.0,
    )
    # Set up the daily context
    ctx = PatientDailyContext(
        hospital_id="hospital_001",
        department=Department.WARD,
        hygiene_level=0.8,
        isolation_effectiveness=0.7,
        diagnostic_speed=0.5,
        is_isolated=False,
        regimen=AntibioticRegimen(on=True, abx_class="beta_lactam", dose_level="std"),
    )
    patient.update_context(ctx)
    return patient


@pytest.fixture
def susceptible_patient() -> Patient:
    """Create a patient in SUSCEPTIBLE state."""
    return Patient(
        patient_id="test_patient_002",
        state=HealthState.SUSCEPTIBLE,
    )


@pytest.fixture
def micro_simulator() -> MicroSimulator:
    """Create a micro simulator with default config."""
    sim = MicroSimulator(config=SimulationConfig())
    yield sim
    sim.close()


# =============================================================================
# Request Structure Tests
# =============================================================================


class TestMicroRequestStructure:
    """Tests that Patient.make_micro_request() produces valid request dicts."""

    def test_carrier_creates_request(self, carrier_patient: Patient):
        """A carrier patient should produce a micro request."""
        request = carrier_patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)
        assert request is not None

    def test_susceptible_returns_none(self, susceptible_patient: Patient):
        """A susceptible patient should NOT produce a micro request."""
        request = susceptible_patient.make_micro_request(
            run_id="run_001", day=1, dt_days=1, seed=42
        )
        assert request is None

    def test_request_has_required_schema_fields(self, carrier_patient: Patient):
        """Request must contain all schema-required top-level fields."""
        request = carrier_patient.make_micro_request(run_id="run_001", day=5, dt_days=1, seed=123)

        required_fields = [
            "schema_version",
            "run_id",
            "episode_id",
            "patient_id",
            "t_day",
            "dt_days",
            "setting",
            "abx",
            "adherence",
            "host",
            "initial_state",
            "seed",
        ]
        for field in required_fields:
            assert field in request, f"Missing required field: {field}"

    def test_request_abx_structure(self, carrier_patient: Patient):
        """ABX sub-dict must have required keys."""
        request = carrier_patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)

        abx = request["abx"]
        assert "on" in abx
        assert "class" in abx
        assert "dose_level" in abx
        assert isinstance(abx["on"], bool)
        assert isinstance(abx["class"], str)
        assert isinstance(abx["dose_level"], str)

    def test_request_host_structure(self, carrier_patient: Patient):
        """Host sub-dict must have required keys."""
        request = carrier_patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)

        host = request["host"]
        required_host_fields = [
            "age_years",
            "immune_strength",
            "immune_status",
            "vulnerability",
            "history_flags",
        ]
        for field in required_host_fields:
            assert field in host, f"Missing host field: {field}"

    def test_request_initial_state_structure(self, carrier_patient: Patient):
        """Initial state sub-dict must have required keys."""
        request = carrier_patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)

        initial = request["initial_state"]
        assert "resistant_fraction" in initial
        assert "dominant_genotype" in initial

    def test_request_types_are_serializable(self, carrier_patient: Patient):
        """All values in request should be JSON-serializable types."""
        import json

        request = carrier_patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)

        # This should not raise
        json_str = json.dumps(request)
        assert json_str is not None

        # Round-trip should preserve structure
        parsed = json.loads(json_str)
        assert parsed == request


# =============================================================================
# Response Structure Tests
# =============================================================================


class TestMicroResponseStructure:
    """Tests that MicroSimulator produces responses compatible with Patient."""

    def test_response_has_required_fields(
        self, carrier_patient: Patient, micro_simulator: MicroSimulator
    ):
        """Response must have all fields that Patient.apply_micro_response expects."""
        request = carrier_patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)
        response = micro_simulator.process_request(request)

        # Required top-level fields
        assert "episode_id" in response
        assert "patient_id" in response
        assert "updated_state" in response
        assert "derived_effects" in response

    def test_response_updated_state_fields(
        self, carrier_patient: Patient, micro_simulator: MicroSimulator
    ):
        """updated_state must contain required fields."""
        request = carrier_patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)
        response = micro_simulator.process_request(request)

        updated_state = response["updated_state"]
        assert "resistant_fraction" in updated_state
        assert "dominant_genotype" in updated_state
        assert "dominant_strain_name" in updated_state

        # Check types
        assert isinstance(updated_state["resistant_fraction"], (int, float))
        assert isinstance(updated_state["dominant_genotype"], str)
        assert isinstance(updated_state["dominant_strain_name"], str)

    def test_response_derived_effects_fields(
        self, carrier_patient: Patient, micro_simulator: MicroSimulator
    ):
        """derived_effects must contain all required fields."""
        request = carrier_patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)
        response = micro_simulator.process_request(request)

        derived = response["derived_effects"]
        required_derived_fields = [
            "relative_transmissibility",
            "p_clearance",
            "severity_modifier",
            "lethality_modifier",
        ]
        for field in required_derived_fields:
            assert field in derived, f"Missing derived_effects field: {field}"
            assert isinstance(derived[field], (int, float)), f"{field} must be numeric"

    def test_response_episode_id_matches_request(
        self, carrier_patient: Patient, micro_simulator: MicroSimulator
    ):
        """Response episode_id must match request episode_id."""
        request = carrier_patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)
        response = micro_simulator.process_request(request)

        assert response["episode_id"] == request["episode_id"]
        assert response["patient_id"] == request["patient_id"]


# =============================================================================
# End-to-End Integration Tests
# =============================================================================


class TestEndToEndIntegration:
    """Tests the full request -> process -> apply cycle."""

    def test_full_cycle_without_errors(
        self, carrier_patient: Patient, micro_simulator: MicroSimulator
    ):
        """Complete request/response cycle should work without exceptions."""
        request = carrier_patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)
        response = micro_simulator.process_request(request)

        # This should not raise
        carrier_patient.apply_micro_response(response)

    def test_patient_state_updated_after_apply(
        self, carrier_patient: Patient, micro_simulator: MicroSimulator
    ):
        """Patient attributes should be updated after apply_micro_response."""

        request = carrier_patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)
        response = micro_simulator.process_request(request)
        carrier_patient.apply_micro_response(response)

        # Values should be updated from response
        assert carrier_patient.resistant_fraction == response["updated_state"]["resistant_fraction"]
        assert carrier_patient.dominant_genotype == response["updated_state"]["dominant_genotype"]
        assert (
            carrier_patient.dominant_strain_name
            == response["updated_state"]["dominant_strain_name"]
        )
        assert (
            carrier_patient.relative_transmissibility
            == response["derived_effects"]["relative_transmissibility"]
        )
        assert carrier_patient.p_clearance == response["derived_effects"]["p_clearance"]

    def test_multi_day_simulation(self, carrier_patient: Patient, micro_simulator: MicroSimulator):
        """Simulate multiple days to ensure state persistence works."""
        for day in range(1, 6):
            request = carrier_patient.make_micro_request(
                run_id="run_001", day=day, dt_days=1, seed=42 + day
            )
            response = micro_simulator.process_request(request)
            carrier_patient.apply_micro_response(response)

            # Patient should still be valid
            assert carrier_patient.state == HealthState.CARRIER
            assert carrier_patient.episode_id == "episode_001"

    def test_response_rejects_mismatched_episode_id(self, carrier_patient: Patient):
        """Patient should reject responses with wrong episode_id."""
        fake_response = {
            "episode_id": "wrong_episode",
            "patient_id": carrier_patient.patient_id,
            "updated_state": {"resistant_fraction": 0.5, "dominant_genotype": "R1"},
            "derived_effects": {
                "relative_transmissibility": 1.0,
                "p_clearance": 0.02,
                "severity_modifier": 1.0,
                "lethality_modifier": 1.0,
            },
        }

        with pytest.raises(ValueError, match="episode_id"):
            carrier_patient.apply_micro_response(fake_response)

    def test_initial_population_uses_dominant_genotype_hint(self):
        """Fresh episodes should seed resistant strains that match the inherited genotype label."""
        population = StrainPopulation.create_initial(
            resistant_fraction=0.6,
            dominant_genotype="R3",
            rng=np.random.default_rng(42),
        )

        _, genotype, _ = get_dominant_strain(population)

        assert genotype == "R3"

    def test_initial_population_assigns_unique_funky_names(self):
        population = StrainPopulation.create_initial(
            resistant_fraction=0.4,
            dominant_genotype="R2",
            rng=np.random.default_rng(7),
        )

        assert len(population.strain_names) == population.n_strains
        assert len(set(population.strain_names)) == population.n_strains
        assert all("_" in name for name in population.strain_names)
        assert all(
            name.split("_", maxsplit=1)[0] != name.split("_", maxsplit=1)[1]
            for name in population.strain_names
        )

    def test_founder_pool_is_deterministic(self):
        config = SimulationConfig(
            founder_pool_size=8,
            founder_pool_seed=11,
            founder_pool_gene_noise_std=0.015,
        )

        pool_a = build_founder_pool(config)
        pool_b = build_founder_pool(config)

        assert [founder.founder_id for founder in pool_a] == [
            founder.founder_id for founder in pool_b
        ]
        assert [founder.founder_name for founder in pool_a] == [
            founder.founder_name for founder in pool_b
        ]
        assert [founder.genotype for founder in pool_a] == [founder.genotype for founder in pool_b]
        for founder_a, founder_b in zip(pool_a, pool_b, strict=False):
            assert np.allclose(founder_a.genome, founder_b.genome)

    def test_initial_population_can_draw_from_founder_pool(self):
        config = SimulationConfig(founder_pool_size=12, founder_pool_seed=1)
        population = StrainPopulation.create_initial(
            resistant_fraction=0.4,
            dominant_genotype="R2",
            rng=np.random.default_rng(7),
            founder_pool=build_founder_pool(config),
            strain_namespace="episode_test",
        )

        assert len(population.strain_ids) == population.n_strains
        assert len(set(population.strain_ids)) == population.n_strains
        assert len(population.parent_ids) == population.n_strains
        assert len(population.founder_ids) == population.n_strains
        assert any(founder_id.startswith("founder_") for founder_id in population.founder_ids)
        assert all(
            strain_id.startswith("episode_test:strain_") for strain_id in population.strain_ids
        )

    def test_episode_state_preserves_strain_names(
        self, carrier_patient: Patient, micro_simulator: MicroSimulator
    ):
        request = carrier_patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)
        response = micro_simulator.process_request(request)
        carrier_patient.apply_micro_response(response)

        state = micro_simulator.get_episode_state("episode_001")
        assert state is not None
        assert len(state.population.strain_names) == state.population.n_strains
        assert carrier_patient.dominant_strain_name in state.population.strain_names
        assert len(state.population.strain_ids) == state.population.n_strains
        assert len(state.population.parent_ids) == state.population.n_strains
        assert len(state.population.donor_ids) == state.population.n_strains
        assert len(state.population.founder_ids) == state.population.n_strains
        assert len(set(state.population.strain_ids)) == state.population.n_strains


# =============================================================================
# Value Range Tests
# =============================================================================


class TestValueRanges:
    """Tests that micro response values are within expected ranges."""

    def test_resistant_fraction_range(
        self, carrier_patient: Patient, micro_simulator: MicroSimulator
    ):
        """resistant_fraction should be in [0, 1]."""
        request = carrier_patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)
        response = micro_simulator.process_request(request)

        rf = response["updated_state"]["resistant_fraction"]
        assert 0.0 <= rf <= 1.0, f"resistant_fraction {rf} out of range [0, 1]"

    def test_p_clearance_range(self, carrier_patient: Patient, micro_simulator: MicroSimulator):
        """p_clearance should be in [0, 1]."""
        request = carrier_patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)
        response = micro_simulator.process_request(request)

        pc = response["derived_effects"]["p_clearance"]
        assert 0.0 <= pc <= 1.0, f"p_clearance {pc} out of range [0, 1]"

    def test_transmissibility_non_negative(
        self, carrier_patient: Patient, micro_simulator: MicroSimulator
    ):
        """relative_transmissibility should be >= 0."""
        request = carrier_patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)
        response = micro_simulator.process_request(request)

        rt = response["derived_effects"]["relative_transmissibility"]
        assert rt >= 0.0, f"relative_transmissibility {rt} should be non-negative"

    def test_modifiers_non_negative(
        self, carrier_patient: Patient, micro_simulator: MicroSimulator
    ):
        """severity_modifier and lethality_modifier should be >= 0."""
        request = carrier_patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)
        response = micro_simulator.process_request(request)

        sm = response["derived_effects"]["severity_modifier"]
        lm = response["derived_effects"]["lethality_modifier"]
        assert sm >= 0.0, f"severity_modifier {sm} should be non-negative"
        assert lm >= 0.0, f"lethality_modifier {lm} should be non-negative"


# =============================================================================
# ABX Scenario Tests
# =============================================================================


class TestABXScenarios:
    """Tests different antibiotic treatment scenarios."""

    @pytest.mark.parametrize(
        "abx_class",
        [
            "none",
            "beta_lactam",
            "fluoroquinolone",
            "aminoglycoside",
            "macrolide",
            "tetracycline",
            "glycopeptide",
        ],
    )
    def test_all_abx_classes_work(self, abx_class: str, micro_simulator: MicroSimulator):
        """All defined ABX classes should be processable."""
        patient = Patient(
            patient_id="test_patient",
            state=HealthState.CARRIER,
            episode_id="episode_001",
        )
        ctx = PatientDailyContext(
            hospital_id="hospital_001",
            department=Department.WARD,
            hygiene_level=0.8,
            isolation_effectiveness=0.7,
            diagnostic_speed=0.5,
            is_isolated=False,
            regimen=AntibioticRegimen(
                on=(abx_class != "none"), abx_class=abx_class, dose_level="std"
            ),
        )
        patient.update_context(ctx)

        request = patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)
        response = micro_simulator.process_request(request)

        # Should complete without error
        patient.apply_micro_response(response)

    @pytest.mark.parametrize("dose_level", ["low", "std", "high"])
    def test_all_dose_levels_work(self, dose_level: str, micro_simulator: MicroSimulator):
        """All dose levels should be processable."""
        patient = Patient(
            patient_id="test_patient",
            state=HealthState.CARRIER,
            episode_id="episode_001",
        )
        ctx = PatientDailyContext(
            hospital_id="hospital_001",
            department=Department.WARD,
            hygiene_level=0.8,
            isolation_effectiveness=0.7,
            diagnostic_speed=0.5,
            is_isolated=False,
            regimen=AntibioticRegimen(on=True, abx_class="beta_lactam", dose_level=dose_level),
        )
        patient.update_context(ctx)

        request = patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)
        response = micro_simulator.process_request(request)

        patient.apply_micro_response(response)


# =============================================================================
# Department Scenario Tests
# =============================================================================


class TestDepartmentScenarios:
    """Tests different department settings."""

    @pytest.mark.parametrize("department", [Department.WARD, Department.ICU])
    def test_all_departments_work(self, department: Department, micro_simulator: MicroSimulator):
        """All department types should be processable."""
        patient = Patient(
            patient_id="test_patient",
            state=HealthState.CARRIER,
            episode_id="episode_001",
        )
        ctx = PatientDailyContext(
            hospital_id="hospital_001",
            department=department,
            hygiene_level=0.8,
            isolation_effectiveness=0.7,
            diagnostic_speed=0.5,
            is_isolated=False,
            regimen=AntibioticRegimen(),
        )
        patient.update_context(ctx)

        request = patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)
        response = micro_simulator.process_request(request)

        patient.apply_micro_response(response)


# =============================================================================
# Immune Status Tests
# =============================================================================


class TestImmuneStatusScenarios:
    """Tests different immune status settings."""

    @pytest.mark.parametrize("immune_status", ["normal", "suppressed"])
    def test_immune_statuses_work(self, immune_status: str, micro_simulator: MicroSimulator):
        """All immune statuses should be processable."""
        patient = Patient(
            patient_id="test_patient",
            state=HealthState.CARRIER,
            episode_id="episode_001",
            immune_status=immune_status,
            immune_strength=0.5 if immune_status == "suppressed" else 1.0,
        )
        ctx = PatientDailyContext(
            hospital_id="hospital_001",
            department=Department.WARD,
            hygiene_level=0.8,
            isolation_effectiveness=0.7,
            diagnostic_speed=0.5,
            is_isolated=False,
            regimen=AntibioticRegimen(),
        )
        patient.update_context(ctx)

        request = patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)
        response = micro_simulator.process_request(request)

        patient.apply_micro_response(response)


# =============================================================================
# Batch Processing Tests
# =============================================================================


class TestBatchProcessing:
    """Tests batch processing of multiple patients."""

    def test_batch_process_multiple_patients(self, micro_simulator: MicroSimulator):
        """Batch processing should handle multiple patients correctly."""
        patients = []
        requests = []

        for i in range(5):
            patient = Patient(
                patient_id=f"patient_{i}",
                state=HealthState.CARRIER,
                episode_id=f"episode_{i}",
            )
            ctx = PatientDailyContext(
                hospital_id="hospital_001",
                department=Department.WARD,
                hygiene_level=0.8,
                isolation_effectiveness=0.7,
                diagnostic_speed=0.5,
                is_isolated=False,
                regimen=AntibioticRegimen(),
            )
            patient.update_context(ctx)
            patients.append(patient)

            request = patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42 + i)
            requests.append(request)

        responses = micro_simulator.process_batch(requests, parallel=False)

        assert len(responses) == len(patients)

        for patient, response in zip(patients, responses):
            patient.apply_micro_response(response)
            # Should all complete without error

    def test_parallel_batch_matches_sequential(self):
        """Parallel batch execution must preserve deterministic results."""
        requests = []
        for i in range(12):
            patient = Patient(
                patient_id=f"patient_{i}",
                state=HealthState.CARRIER,
                episode_id=f"episode_{i}",
            )
            ctx = PatientDailyContext(
                hospital_id="hospital_001",
                department=Department.WARD,
                hygiene_level=0.8,
                isolation_effectiveness=0.7,
                diagnostic_speed=0.5,
                is_isolated=False,
                regimen=AntibioticRegimen(on=True, abx_class="beta_lactam", dose_level="std"),
            )
            patient.update_context(ctx)
            requests.append(
                patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42 + i)
            )

        seq = MicroSimulator(config=SimulationConfig(), n_workers=2)
        par = MicroSimulator(config=SimulationConfig(), n_workers=2)
        try:
            seq_responses = seq.process_batch(requests, parallel=False)
            par_responses = par.process_batch(requests, parallel=True)
        finally:
            seq.close()
            par.close()

        assert par_responses == seq_responses

    def test_duplicate_episode_ids_in_batch_raise(self, micro_simulator: MicroSimulator):
        """The same episode cannot be advanced twice in one parallel batch."""
        patient = Patient(
            patient_id="patient_dup",
            state=HealthState.CARRIER,
            episode_id="episode_dup",
        )
        ctx = PatientDailyContext(
            hospital_id="hospital_001",
            department=Department.WARD,
            hygiene_level=0.8,
            isolation_effectiveness=0.7,
            diagnostic_speed=0.5,
            is_isolated=False,
            regimen=AntibioticRegimen(),
        )
        patient.update_context(ctx)
        request = patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)

        with pytest.raises(ValueError, match="Duplicate episode_id"):
            micro_simulator.process_batch([request, dict(request)], parallel=True)


# =============================================================================
# Edge Case Tests
# =============================================================================


class TestEdgeCases:
    """Tests edge cases and boundary conditions."""

    def test_zero_adherence(self, micro_simulator: MicroSimulator):
        """Patient with zero adherence should be processable."""
        patient = Patient(
            patient_id="test_patient",
            state=HealthState.CARRIER,
            episode_id="episode_001",
            compliance=0.0,  # Zero compliance -> zero adherence
        )
        ctx = PatientDailyContext(
            hospital_id="hospital_001",
            department=Department.WARD,
            hygiene_level=0.8,
            isolation_effectiveness=0.7,
            diagnostic_speed=0.5,
            is_isolated=False,
            regimen=AntibioticRegimen(on=True, abx_class="beta_lactam"),
        )
        patient.update_context(ctx)

        request = patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)
        response = micro_simulator.process_request(request)
        patient.apply_micro_response(response)

    def test_high_vulnerability_patient(self, micro_simulator: MicroSimulator):
        """Patient with very high vulnerability should be processable."""
        patient = Patient(
            patient_id="test_patient",
            state=HealthState.CARRIER,
            episode_id="episode_001",
            vulnerability=3.0,
            immune_strength=0.3,
        )
        ctx = PatientDailyContext(
            hospital_id="hospital_001",
            department=Department.ICU,
            hygiene_level=0.9,
            isolation_effectiveness=0.9,
            diagnostic_speed=0.8,
            is_isolated=True,
            regimen=AntibioticRegimen(on=True, abx_class="glycopeptide", dose_level="high"),
        )
        patient.update_context(ctx)

        request = patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)
        response = micro_simulator.process_request(request)
        patient.apply_micro_response(response)

    def test_patient_with_history_flags(self, micro_simulator: MicroSimulator):
        """Patient with history flags should be processable."""
        patient = Patient(
            patient_id="test_patient",
            state=HealthState.CARRIER,
            episode_id="episode_001",
            history_flags={"prior_abx", "diabetes", "recent_surgery", "immunosuppression"},
        )
        ctx = PatientDailyContext(
            hospital_id="hospital_001",
            department=Department.WARD,
            hygiene_level=0.8,
            isolation_effectiveness=0.7,
            diagnostic_speed=0.5,
            is_isolated=False,
            regimen=AntibioticRegimen(),
        )
        patient.update_context(ctx)

        request = patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)

        # History flags should be present and sorted
        assert "history_flags" in request["host"]
        assert request["host"]["history_flags"] == sorted(patient.history_flags)

        response = micro_simulator.process_request(request)
        patient.apply_micro_response(response)

    def test_missing_context_raises_error(self):
        """Patient without context should raise error on make_micro_request."""
        patient = Patient(
            patient_id="test_patient",
            state=HealthState.CARRIER,
            episode_id="episode_001",
        )
        # Don't call update_context()

        with pytest.raises(RuntimeError, match="PatientDailyContext not set"):
            patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)


# =============================================================================
# Reproducibility Tests
# =============================================================================


class TestReproducibility:
    """Tests that simulations are reproducible with same seed."""

    def test_same_seed_same_result(self, micro_simulator: MicroSimulator):
        """Same seed should produce identical results."""

        def run_simulation(seed: int) -> dict:
            patient = Patient(
                patient_id="test_patient",
                state=HealthState.CARRIER,
                episode_id="episode_001",
            )
            ctx = PatientDailyContext(
                hospital_id="hospital_001",
                department=Department.WARD,
                hygiene_level=0.8,
                isolation_effectiveness=0.7,
                diagnostic_speed=0.5,
                is_isolated=False,
                regimen=AntibioticRegimen(on=True, abx_class="beta_lactam"),
            )
            patient.update_context(ctx)

            request = patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=seed)
            # Use a fresh simulator to avoid state pollution
            sim = MicroSimulator(config=SimulationConfig())
            return sim.process_request(request)

        result1 = run_simulation(seed=12345)
        result2 = run_simulation(seed=12345)

        # Results should be identical
        assert (
            result1["updated_state"]["resistant_fraction"]
            == result2["updated_state"]["resistant_fraction"]
        )
        assert (
            result1["updated_state"]["dominant_genotype"]
            == result2["updated_state"]["dominant_genotype"]
        )
        assert (
            result1["derived_effects"]["p_clearance"] == result2["derived_effects"]["p_clearance"]
        )

    def test_different_seed_different_result(self, micro_simulator: MicroSimulator):
        """Different seeds should generally produce different results."""

        def run_simulation(seed: int) -> dict:
            patient = Patient(
                patient_id="test_patient",
                state=HealthState.CARRIER,
                episode_id="episode_001",
            )
            ctx = PatientDailyContext(
                hospital_id="hospital_001",
                department=Department.WARD,
                hygiene_level=0.8,
                isolation_effectiveness=0.7,
                diagnostic_speed=0.5,
                is_isolated=False,
                regimen=AntibioticRegimen(on=True, abx_class="beta_lactam"),
            )
            patient.update_context(ctx)

            request = patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=seed)
            sim = MicroSimulator(config=SimulationConfig())
            return sim.process_request(request)

        results = [run_simulation(seed=i) for i in range(10)]

        # At least some results should differ (stochastic simulation)
        clearances = [r["derived_effects"]["p_clearance"] for r in results]
        # With 10 different seeds, we expect some variation
        assert len(set(clearances)) > 1, "All seeds produced identical p_clearance values"


# =============================================================================
# Lifecycle / Stochasticity Tests
# =============================================================================


class TestLifecycleDynamics:
    """Tests lineage aging, genetic buffering, and stochastic extinction."""

    def test_episode_state_tracks_lifecycle_fields(self, carrier_patient: Patient):
        sim = MicroSimulator(config=SimulationConfig())

        request = carrier_patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)
        sim.process_request(request)

        state = sim.get_episode_state("episode_001")
        assert state is not None
        assert len(state.population.lineage_ages) == state.population.n_strains
        assert len(state.population.damage_loads) == state.population.n_strains
        assert np.all(state.population.lineage_ages >= 0.0)
        assert np.all(state.population.damage_loads >= 0.0)

    def test_mutations_record_parent_strain_ids(self, carrier_patient: Patient):
        sim = MicroSimulator(
            config=SimulationConfig(
                founder_pool_size=10,
                base_mutation_rate=1.0,
                mutation_std=0.01,
                base_hgt_rate=0.0,
                max_strains=80,
            )
        )

        request = carrier_patient.make_micro_request(run_id="run_001", day=1, dt_days=1, seed=42)
        sim.process_request(request)

        state = sim.get_episode_state("episode_001")
        assert state is not None
        assert any(parent_id for parent_id in state.population.parent_ids)

    def test_repair_and_persistence_reduce_turnover_costs(self):
        fragile = create_wild_type_genome()
        fragile[GeneIndex.DNA_REPAIR] = 0.0
        fragile[GeneIndex.DORMANCY_PROPENSITY] = 0.0
        fragile[GeneIndex.STRESS_RESPONSE] = 0.0
        fragile[GeneIndex.DAMAGE_TOLERANCE] = 0.0

        buffered = create_wild_type_genome()
        buffered[GeneIndex.DNA_REPAIR] = 1.0
        buffered[GeneIndex.DORMANCY_PROPENSITY] = 1.0
        buffered[GeneIndex.STRESS_RESPONSE] = 1.0
        buffered[GeneIndex.DAMAGE_TOLERANCE] = 1.0

        population = StrainPopulation(
            genomes=np.vstack([fragile, buffered]).astype(np.float32),
            populations=np.array([5e4, 5e4], dtype=np.float64),
        )

        final_pop, _ = simulate_day(
            population=population,
            abx_class="beta_lactam",
            dose_level="std",
            adherence=1.0,
            immune_strength=0.8,
            immune_status="normal",
            config=SimulationConfig(
                base_mutation_rate=0.0,
                base_hgt_rate=0.0,
                stochastic_threshold=0.0,
                stochastic_noise_scale=0.0,
            ),
            seed=7,
        )

        assert final_pop.damage_loads[1] < final_pop.damage_loads[0]
        assert final_pop.populations[1] > final_pop.populations[0]

    def test_gene_combinations_outperform_single_specialists(self):
        repair_only = create_wild_type_genome()
        repair_only[GeneIndex.DNA_REPAIR] = 1.0
        repair_only[GeneIndex.DORMANCY_PROPENSITY] = 0.0
        repair_only[GeneIndex.STRESS_RESPONSE] = 0.0
        repair_only[GeneIndex.DAMAGE_TOLERANCE] = 0.0

        dormancy_only = create_wild_type_genome()
        dormancy_only[GeneIndex.DNA_REPAIR] = 0.0
        dormancy_only[GeneIndex.DORMANCY_PROPENSITY] = 1.0
        dormancy_only[GeneIndex.STRESS_RESPONSE] = 0.0
        dormancy_only[GeneIndex.DAMAGE_TOLERANCE] = 0.0

        combo = create_wild_type_genome()
        combo[GeneIndex.DNA_REPAIR] = 1.0
        combo[GeneIndex.DORMANCY_PROPENSITY] = 1.0
        combo[GeneIndex.STRESS_RESPONSE] = 1.0
        combo[GeneIndex.DAMAGE_TOLERANCE] = 1.0

        population = StrainPopulation(
            genomes=np.vstack([repair_only, dormancy_only, combo]).astype(np.float32),
            populations=np.array([4e4, 4e4, 4e4], dtype=np.float64),
        )

        final_pop, _ = simulate_day(
            population=population,
            abx_class="beta_lactam",
            dose_level="std",
            adherence=1.0,
            immune_strength=0.8,
            immune_status="normal",
            config=SimulationConfig(
                base_mutation_rate=0.0,
                base_hgt_rate=0.0,
                stochastic_threshold=0.0,
                stochastic_noise_scale=0.0,
            ),
            seed=17,
        )

        assert final_pop.populations[2] > final_pop.populations[0]
        assert final_pop.populations[2] > final_pop.populations[1]

    def test_small_population_can_die_out_stochastically(self):
        population = StrainPopulation(
            genomes=create_wild_type_genome().reshape(1, -1),
            populations=np.array([5.0], dtype=np.float64),
        )

        final_pop, _ = simulate_day(
            population=population,
            abx_class="glycopeptide",
            dose_level="high",
            adherence=1.0,
            immune_strength=2.0,
            immune_status="normal",
            config=SimulationConfig(base_mutation_rate=0.0, base_hgt_rate=0.0),
            seed=123,
        )

        assert final_pop.total_population == 0.0
