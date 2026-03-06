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
    base_diagnostic_speed: float = 0.5  # ≥0; higher ⇒ faster detection

    # Transmission
    base_transmission_rate: float = 0.05  # per-carrier-per-susceptible daily beta

    # Antibiotic policy (daily probability that regimen.on=True)
    icu_abx_probability: float = 0.60
    ward_abx_probability: float = 0.15

    # Detection / isolation of carriers
    carrier_isolation_probability: float = 0.30  # daily detection chance for non-isolated carriers
