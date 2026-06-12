"""Shared data models for the micro simulation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class FounderStrain:
    """Reusable founder strain used to initialize patient episodes.

    Attributes:
        founder_id: Stable identifier for the founder record.
        founder_name: Human-readable strain name.
        genotype: Resistance class label, such as ``S``, ``R1``, ``R2``, or ``R3``.
        genome: Founder genome vector with one value per micro gene.
    """

    founder_id: str
    founder_name: str
    genotype: str
    genome: np.ndarray
