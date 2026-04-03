from __future__ import annotations

import mesa

from mss.domain import Patient


class _MesaModel(mesa.Model):
    """Minimal Mesa model stub for agent construction."""


class PatientAgent(mesa.Agent):
    """Mesa-compatible wrapper around the exchange-layer Patient."""

    def __init__(self, patient: Patient, model: mesa.Model) -> None:
        super().__init__(model)
        self.patient = patient

    def step(self) -> None:
        return
