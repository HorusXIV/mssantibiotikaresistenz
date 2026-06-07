from __future__ import annotations

import mesa

from mss.domain import Patient


class _MesaModel(mesa.Model):
    """Minimal Mesa model wrapper."""


class PatientAgent(mesa.Agent):
    """Mesa-compatible wrapper around the exchange-layer Patient."""

    def __init__(self, patient: Patient, model: mesa.Model) -> None:
        super().__init__(model)
        self.patient = patient

    def step(self) -> None:
        # No per-agent logic; the macro simulator drives patient state.
        return
