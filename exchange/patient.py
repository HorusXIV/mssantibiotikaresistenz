from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, Optional, Set


# -------------------------
# Small enums for clarity
# -------------------------

class HealthState(str, Enum):
    SUSCEPTIBLE = "S"
    CARRIER = "C"


class Department(str, Enum):
    WARD = "ward"
    ICU = "icu"


class TreatmentPhase(str, Enum):
    NONE = "none"
    ON_ABX = "on_abx"
    POST_ABX = "post_abx"


# -------------------------
# Context snapshots (macro -> patient)
# -------------------------

@dataclass
class AntibioticRegimen:
    on: bool = False
    abx_class: str = "none"          # e.g. "beta_lactam", "fluoroquinolone"
    dose_level: str = "std"          # "low"|"std"|"high" (or numeric)
    # If you later want duration tracking, keep it in macro or add days_on_abx here.


@dataclass
class PatientDailyContext:
    """
    A compact snapshot of the environment the patient is in today.
    This comes from the macro/hospital layer.
    """
    hospital_id: str
    department: Department

    # Location-level controls (macro-owned)
    hygiene_level: float             # 0..1 (1 = best hygiene). Macro decides semantics.
    isolation_effectiveness: float   # 0..1 (1 = perfect isolation)
    diagnostic_speed: float          # e.g. tests/day or inverse delay; macro decides semantics.

    # Patient-specific in macro
    is_isolated: bool
    regimen: AntibioticRegimen


# -------------------------
# Patient class (the interface object)
# -------------------------

@dataclass
class Patient:
    # Identity and macro state
    patient_id: str
    state: HealthState = HealthState.SUSCEPTIBLE
    episode_id: Optional[str] = None

    # Where they currently are (macro will overwrite via update_context)
    hospital_id: Optional[str] = None
    department: Department = Department.WARD
    is_isolated: bool = False

    # Intrinsic "base stats" (keep these stable during an episode)
    age_years: int = 50
    compliance: float = 0.8          # 0..1 (affects adherence, isolation compliance, etc.)
    vulnerability: float = 1.0       # multiplier >= ~0.5 .. 3.0 (your choice)
    immune_strength: float = 1.0     # multiplier; 1.0 = normal
    immune_status: str = "normal"    # "normal"|"suppressed" (or more categories)
    sociability: float = 1.0         # contact-rate multiplier (macro uses)

    # Treatment info (patient-specific)
    treatment_phase: TreatmentPhase = TreatmentPhase.NONE
    adherence: float = 1.0           # 0..1, can be derived from compliance each day
    regimen: AntibioticRegimen = field(default_factory=AntibioticRegimen)

    # Medical history (keep it simple: flags)
    history_flags: Set[str] = field(default_factory=set)
    # Examples: {"prior_abx", "diabetes", "recent_surgery", "immunosuppression"}

    # --- Micro outputs stored on patient (updated daily if carrier) ---
    resistant_fraction: float = 0.0              # 0..1
    dominant_genotype: str = "S"                 # "S", "R1", "R2", ...
    relative_transmissibility: float = 1.0       # multiplies macro beta
    p_clearance: float = 0.02                    # daily probability C->S
    severity_modifier: float = 1.0              # strain-driven (optional for macro)
    lethality_modifier: float = 1.0             # strain-driven (optional for macro)

    # Store last daily context for reproducibility/debugging
    _ctx: Optional[PatientDailyContext] = field(default=None, repr=False)

    # -------------------------
    # Macro -> Patient
    # -------------------------
    def update_context(self, ctx: PatientDailyContext) -> None:
        """
        Called by macro once per day before micro requests are built.
        """
        self._ctx = ctx
        self.hospital_id = ctx.hospital_id
        self.department = ctx.department
        self.is_isolated = ctx.is_isolated
        self.regimen = ctx.regimen

        # Derive adherence from compliance (simple mapping; tweak as you like)
        # ICU might have better adherence due to supervision; you can encode that here too.
        base = self.compliance
        if self.department == Department.ICU:
            base = min(1.0, base + 0.1)
        self.adherence = max(0.0, min(1.0, base))

        # Treatment phase (for micro; keep very simple)
        if self.regimen.on:
            self.treatment_phase = TreatmentPhase.ON_ABX
        elif self.treatment_phase == TreatmentPhase.ON_ABX and not self.regimen.on:
            self.treatment_phase = TreatmentPhase.POST_ABX
        elif not self.regimen.on:
            self.treatment_phase = TreatmentPhase.NONE

    # -------------------------
    # Patient -> Micro request (daily if carrier)
    # -------------------------
    def make_micro_request(self, run_id: str, day: int, dt_days: int, seed: int) -> Optional[Dict[str, Any]]:
        """
        Returns a dict that your MicroAdapter writes to JSONL.
        Only returns a request if patient is a carrier.
        """
        if self.state != HealthState.CARRIER or self.episode_id is None:
            return None
        if self._ctx is None:
            raise RuntimeError("PatientDailyContext not set. Call update_context() first.")

        return {
            "schema_version": "1.0",
            "run_id": run_id,
            "episode_id": self.episode_id,
            "patient_id": self.patient_id,
            "t_day": day,
            "dt_days": dt_days,

            # Environment / selection inputs
            "setting": self.department.value,
            "abx": {
                "on": self.regimen.on,
                "class": self.regimen.abx_class,
                "dose_level": self.regimen.dose_level,
            },
            "adherence": self.adherence,

            # Host factors (micro environment)
            "host": {
                "age_years": self.age_years,
                "immune_strength": self.immune_strength,
                "immune_status": self.immune_status,
                "vulnerability": self.vulnerability,
                "history_flags": sorted(list(self.history_flags)),
            },

            # State handed over to micro
            "initial_state": {
                "resistant_fraction": self.resistant_fraction,
                "dominant_genotype": self.dominant_genotype,
            },

            # Reproducibility
            "seed": seed,
        }

    # -------------------------
    # Micro -> Patient response
    # -------------------------
    def apply_micro_response(self, resp: Dict[str, Any]) -> None:
        """
        Applies a micro response (one episode/patient/day).
        The macro layer remains the "owner" of S/C transitions,
        but micro can supply p_clearance and multipliers.
        """
        # Minimal safety checks
        if resp.get("episode_id") != self.episode_id:
            raise ValueError("Micro response episode_id does not match patient episode_id.")

        updated_state = resp.get("updated_state", {})
        derived = resp.get("derived_effects", {})

        # Within-host state
        if "resistant_fraction" in updated_state:
            self.resistant_fraction = float(updated_state["resistant_fraction"])
        if "dominant_genotype" in updated_state:
            self.dominant_genotype = str(updated_state["dominant_genotype"])

        # Derived effects for macro
        if "relative_transmissibility" in derived:
            self.relative_transmissibility = float(derived["relative_transmissibility"])
        if "p_clearance" in derived:
            self.p_clearance = float(derived["p_clearance"])

        # Optional strain-driven outcomes
        if "severity_modifier" in derived:
            self.severity_modifier = float(derived["severity_modifier"])
        if "lethality_modifier" in derived:
            self.lethality_modifier = float(derived["lethality_modifier"])

    # -------------------------
    # Helper getters for macro
    # -------------------------
    def transmission_multiplier_for_macro(self) -> float:
        """
        Macro can multiply its base beta by this.
        Includes patient sociability and micro strain transmissibility.
        Isolation/hygiene are handled at hospital-level, not here.
        """
        return max(0.0, self.sociability) * max(0.0, self.relative_transmissibility)

    def susceptibility_multiplier_for_macro(self) -> float:
        """
        Used for S->C infection/colonization probability in macro.
        Keep it simple: vulnerability and immune_strength could influence it.
        """
        # Higher immune_strength reduces susceptibility.
        immune_effect = 1.0 / max(0.1, self.immune_strength)
        return max(0.0, self.vulnerability) * immune_effect

    def daily_death_risk_multiplier(self) -> float:
        """
        Only relevant if you model death in macro.
        This is a multiplier on a base death risk.
        """
        # Example: vulnerability and strain lethality
        return max(0.0, self.vulnerability) * max(0.0, self.lethality_modifier)

    def should_clear_today(self, rng) -> bool:
        """
        Macro calls this to decide C->S.
        """
        return (self.state == HealthState.CARRIER) and (rng.random() < self.p_clearance)

    def clear_carriage(self) -> None:
        """
        Macro-owned transition: C -> S
        """
        self.state = HealthState.SUSCEPTIBLE
        self.episode_id = None
        self.resistant_fraction = 0.0
        self.dominant_genotype = "S"
        self.relative_transmissibility = 1.0
        self.p_clearance = 0.02
        self.severity_modifier = 1.0
        self.lethality_modifier = 1.0