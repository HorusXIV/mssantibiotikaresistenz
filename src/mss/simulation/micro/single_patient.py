"""Data models for the single-patient micro runner."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class ABXPeriod:
    """Antibiotic exposure interval for a single-patient micro run.

    Attributes:
        start_day: First simulated day on which the antibiotic is active.
        end_day: Last simulated day on which the antibiotic is active.
        abx_class: Antibiotic class passed to the micro engine.
        dose_level: Dose label passed to the micro engine.
        adherence: Effective adherence multiplier in ``[0, 1]``.
    """

    start_day: int
    end_day: int
    abx_class: str
    dose_level: str
    adherence: float


@dataclass
class EpisodeConfig:
    """Initial conditions for a single-patient micro episode.

    Attributes:
        n_days: Number of macro days to simulate.
        seed: Base random seed for deterministic replay.
        resistant_fraction: Initial population fraction in resistant strains.
        dominant_genotype: Initial dominant genotype class.
        initial_population: Initial within-host bacterial population size.
        immune_strength: Host immune-strength multiplier.
        abx_schedule: Antibiotic exposure periods to apply during the episode.
        allow_spontaneous_clearance: Whether to roll for spontaneous clearance each day.
    """

    n_days: int
    seed: int
    resistant_fraction: float
    dominant_genotype: str
    initial_population: float
    immune_strength: float
    abx_schedule: list[ABXPeriod]
    allow_spontaneous_clearance: bool = True


@dataclass
class DayRecord:
    """Daily output record from a single-patient micro episode.

    Attributes:
        day: Simulated macro day.
        total_population: Total bacterial population at the end of the day.
        resistant_fraction: Population-weighted resistant fraction.
        n_strains: Number of active strains after consolidation.
        shannon_entropy: Shannon diversity of strain populations.
        p_clearance: Daily clearance probability derived from micro state.
        abx_class: Antibiotic class active on this day, or ``none``.
        mean_genes: Population-weighted mean genome vector.
        mean_damage: Population-weighted mean damage load.
        mean_age: Population-weighted mean lineage age.
        frac_S: Population fraction classified as susceptible.
        frac_R1: Population fraction classified as low resistance.
        frac_R2: Population fraction classified as medium resistance.
        frac_R3: Population fraction classified as high resistance.
        strain_snapshot: Per-strain tuples of name, population, and genotype class.
        cleared: Whether the population has fallen below the clearance floor.
    """

    day: int
    total_population: float
    resistant_fraction: float
    n_strains: int
    shannon_entropy: float
    p_clearance: float
    abx_class: str
    mean_genes: np.ndarray
    mean_damage: float
    mean_age: float
    frac_S: float
    frac_R1: float
    frac_R2: float
    frac_R3: float
    strain_snapshot: list[tuple[str, float, str]]
    cleared: bool = False
