"""
Simulation configuration for the macro (hospital/network) layer.

The macro layer models inter-patient transmission within and across hospitals
using a simple carrier model (states: SUSCEPTIBLE and CARRIER).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class SimulationConfig:
    """Parameters that govern macro-level simulation behaviour.

    All probabilities are daily and lie in [0, 1].
    """

    # Hospital environment
    base_hygiene: float = 0.7  # 0..1; higher ⇒ less transmission
    base_isolation_effectiveness: float = 0.8  # 0..1; how well isolation works
    base_diagnostic_speed: float = 1.0  # >0; multiplier on carrier detection (1.0 = unchanged)

    # Transmission
    base_transmission_rate: float = 0.05  # per-carrier-per-susceptible daily beta
    daily_contact_attempts: float = 80.0  # effective daily patient contacts for hazard model

    # Antibiotic policy (daily probability that regimen.on=True)
    icu_abx_probability: float = 0.60
    ward_abx_probability: float = 0.15

    # Detection / isolation of carriers
    carrier_isolation_probability: float = 0.30  # daily detection chance for non-isolated carriers

    # LOS (length of stay): log-normal distribution drawn at admission
    los_mean_ward: float = 8.0  # mean ward stay in days
    los_mean_icu: float = 14.0  # mean ICU stay in days
    los_sigma: float = 0.6  # log-normal sigma (shape parameter)

    # Discharge
    carrier_extension_days: float = 14.0  # extend stay by this many days when carrier hits target
    discharge_logistic_k: float = 1.0  # steepness of logistic discharge curve
    discharge_logistic_t_half: float = 3.0  # days over target where p=50%

    # Admissions (0 = disabled)
    daily_admission_rate: float = 0.0  # total new admissions per day (Poisson mean)
    community_carrier_fraction: float = 0.03
    replacement_resistant_fraction: float = 0.05
    replacement_dominant_genotype: str = "S"
    max_occupancy_per_hospital: int = 200  # hard cap per hospital

    # Daily mortality (carriers scaled by severity_modifier). Derived from CH
    # in-hospital mortality ~2.5%/admission over ~5.5-day LOS (BFS 2022).
    base_mortality_rate: float = 0.0045

    # Inter-hospital transfers (0.0 = disabled)
    # Per-patient daily probability of being transferred; destination is chosen
    # from all other hospitals, weighted by (free_capacity / distance).
    daily_transfer_rate: float = 0.0

    # Grid
    dept_grid_cols: int = 3
    dept_grid_rows: int = 2
    dept_grid_icu_rows: int = 1
    network_grid_cols: int = 3
    proximity_decay_alpha: float = 0.5
