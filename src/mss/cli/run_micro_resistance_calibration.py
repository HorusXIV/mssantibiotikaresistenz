"""Code-level calibration of within-host resistance persistence (Calibration 6).

The micro layer washes out the resistant fraction within days without antibiotics,
which contradicts the months-long MRSA carriage observed clinically. The dominant
driver is the per-step fitness gap between sensitive and resistant strains, which is
set by *code constants*, not by the scenario YAML:

  (a) ``ResistanceCosts`` defaults in ``micro/genome.py``
  (b) the resistant seed-genome ``GROWTH_BASE`` in
      ``micro/engine.py::_create_seed_genome_for_genotype``

This tool sweeps two interpretable levers over those constants and reports, for each
grid point, how long the within-host resistant fraction survives without antibiotics
(``rf_half_life_days``: days to fall from 0.90 to 0.45). It also runs an antibiotic
counter-check so a candidate is only accepted if resistance still gets *selected*
(rather than the population going extinct) under beta-lactam pressure.

It uses only the micro engine (no macro/mesa), so it is cheap and dependency-light.

Usage:
    uv run python -m mss.cli.run_micro_resistance_calibration \
        --config config/simulation_realistic_micro.yml
    uv run python -m mss.cli.run_micro_resistance_calibration \
        --config config/simulation_realistic_micro.yml \
        --cost-scale 1.0 0.6 0.3 0.15 --growth-gap 0.12 0.08 0.04 0.0

Levers:
    --cost-scale   multiplier on ResistanceCosts (1.0 = current code, <1 = compensated)
    --growth-gap   GROWTH_BASE difference S minus resistant (0.12 = current code, 0 = none)

Outputs (outputs/<timestamp>_MicroResistance/data/):
    micro_resistance_grid.parquet          one row per (cost_scale, growth_gap)
    micro_resistance_trajectories.parquet  daily resistant_fraction per grid point
"""

from __future__ import annotations

import argparse
from dataclasses import replace
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

import mss.simulation.micro.engine as engine_module
import mss.simulation.micro.genome as genome_module
from mss.cli.run_coupled_simulation import PROJECT_ROOT, load_coupled_settings
from mss.simulation.micro.engine import (
    SimulationConfig,
    StrainPopulation,
    simulate_day,
)
from mss.simulation.micro.genome import GeneIndex, compute_resistant_fraction

# Current-code reference values for the two levers (kept here so the patches are
# self-documenting and reversible).
_REFERENCE_COSTS = dict(
    efflux_pumps=0.15,
    target_modification=0.12,
    permeability_reduction=0.08,
)
_SENSITIVE_GROWTH_BASE = 0.80  # create_wild_type_genome() GROWTH_BASE

_ORIGINAL_RESISTANCE_COSTS_INIT = genome_module.ResistanceCosts.__init__
_ORIGINAL_SEED_GENOME = engine_module._create_seed_genome_for_genotype


def _apply_cost_scale(scale: float) -> None:
    """Override ResistanceCosts() defaults so the engine sees scaled costs."""
    eff = _REFERENCE_COSTS["efflux_pumps"] * scale
    tgt = _REFERENCE_COSTS["target_modification"] * scale
    perm = _REFERENCE_COSTS["permeability_reduction"] * scale

    def __init__(
        self,
        efflux_pumps: float = eff,
        target_modification: float = tgt,
        permeability_reduction: float = perm,
    ) -> None:
        self.efflux_pumps = efflux_pumps
        self.target_modification = target_modification
        self.permeability_reduction = permeability_reduction

    genome_module.ResistanceCosts.__init__ = __init__


def _apply_growth_gap(gap: float) -> None:
    """Override resistant seed-genome GROWTH_BASE to S minus ``gap``."""
    target_growth = float(np.clip(_SENSITIVE_GROWTH_BASE - gap, 0.0, 1.0))

    def patched(dominant_genotype: str) -> np.ndarray:
        genome = _ORIGINAL_SEED_GENOME(dominant_genotype)
        if dominant_genotype in ("R1", "R2", "R3"):
            genome[GeneIndex.GROWTH_BASE] = target_growth
        return genome

    engine_module._create_seed_genome_for_genotype = patched


def _restore_levers() -> None:
    genome_module.ResistanceCosts.__init__ = _ORIGINAL_RESISTANCE_COSTS_INIT
    engine_module._create_seed_genome_for_genotype = _ORIGINAL_SEED_GENOME


def _resistant_fraction_trajectory(
    config: SimulationConfig,
    *,
    abx_on: bool,
    abx_class: str,
    dose_level: str,
    adherence: float,
    immune_strength: float,
    resistant_fraction: float,
    dominant_genotype: str,
    initial_population: float,
    n_days: int,
    n_seeds: int,
) -> tuple[np.ndarray, list[int | None]]:
    """Run controlled single-carrier episodes; return mean rf series and clear days."""
    abx = abx_class if abx_on else "none"
    series: list[list[float]] = []
    clear_days: list[int | None] = []
    for seed in range(n_seeds):
        rng = np.random.default_rng(seed)
        population = StrainPopulation.create_initial(
            resistant_fraction=resistant_fraction,
            dominant_genotype=dominant_genotype,
            initial_population=initial_population,
            rng=rng,
            strain_namespace=f"cal6_{seed}",
        )
        rf_series: list[float] = []
        cleared: int | None = None
        for day in range(1, n_days + 1):
            population, _ = simulate_day(
                population=population,
                abx_class=abx,
                dose_level=dose_level,
                adherence=adherence,
                immune_strength=immune_strength,
                config=config,
                seed=seed * 100_000 + day,
            )
            rf_series.append(compute_resistant_fraction(population.genomes, population.populations))
            if cleared is None and population.total_population < config.min_population:
                cleared = day
        series.append(rf_series)
        clear_days.append(cleared)
    return np.asarray(series).mean(axis=0), clear_days


def _half_life_days(mean_rf: np.ndarray, start: float = 0.90, target: float = 0.45) -> float:
    """Days until the mean resistant fraction first drops to/below ``target``.

    Returns ``n_days`` (the run length) if it never drops that far, signalling a
    stable plateau that meets the persistence goal.
    """
    for day_index, value in enumerate(mean_rf, start=1):
        if value <= target:
            return float(day_index)
    return float(len(mean_rf))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calibrate code-level resistance persistence (Calibration 6)."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=PROJECT_ROOT / "config" / "simulation_realistic_micro.yml",
    )
    parser.add_argument(
        "--cost-scale",
        type=float,
        nargs="+",
        default=[1.0, 0.6, 0.3, 0.15],
        help="Multipliers on ResistanceCosts (1.0 = current code).",
    )
    parser.add_argument(
        "--growth-gap",
        type=float,
        nargs="+",
        default=[0.12, 0.08, 0.04, 0.0],
        help="GROWTH_BASE gap S minus resistant (0.12 = current code).",
    )
    parser.add_argument(
        "--death-floor",
        type=float,
        nargs="+",
        default=[0.1, 1.0],
        help=(
            "micro.death_fitness_floor values. Floor in death=death_rate/(fitness+floor); "
            "0.1 = legacy hidden selection channel, higher = weaker. Config-reachable."
        ),
    )
    parser.add_argument(
        "--selection-strength",
        type=float,
        nargs="+",
        default=None,
        help=(
            "Optional micro.selection_strength values to co-sweep. If omitted, the "
            "config's selection_strength is used unchanged."
        ),
    )
    parser.add_argument("--n-days", type=int, default=120)
    parser.add_argument("--n-seeds", type=int, default=8)
    parser.add_argument("--resistant-fraction", type=float, default=0.90)
    parser.add_argument("--dominant-genotype", type=str, default="R2")
    parser.add_argument("--initial-population", type=float, default=1e6)
    parser.add_argument("--immune-strength", type=float, default=0.75)
    parser.add_argument("--adherence", type=float, default=0.70)
    parser.add_argument("--target-half-life-days", type=float, default=90.0)
    parser.add_argument("--counter-check-abx-class", type=str, default="beta_lactam")
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    settings = load_coupled_settings(args.config)
    base_config = settings.micro

    output_dir = args.output_dir or (
        PROJECT_ROOT / "outputs" / (datetime.now().strftime("%Y%m%d_%H%M%S") + "_MicroResistance")
    )
    data_dir = output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print("Calibration 6: within-host resistance persistence (code + config levers)")
    print(f"  Config: {args.config}")
    print(f"  cost_scale grid:     {args.cost_scale}")
    print(f"  growth_gap grid:     {args.growth_gap}")
    print(f"  death_floor grid:    {args.death_floor}")
    print(
        "  selection_strength:  "
        + (
            str(args.selection_strength)
            if args.selection_strength
            else f"[config={base_config.selection_strength:g}]"
        )
    )
    print(
        f"  Episode: rf0={args.resistant_fraction}, {args.dominant_genotype}, "
        f"{args.n_days} days, {args.n_seeds} seeds, no ABX"
    )
    print(f"  Target: rf_half_life_days >= {args.target_half_life_days:g}")
    print("=" * 72)

    selection_values = (
        args.selection_strength if args.selection_strength else [base_config.selection_strength]
    )

    grid_rows: list[dict[str, Any]] = []
    traj_rows: list[dict[str, Any]] = []

    try:
        for cost_scale in args.cost_scale:
            for growth_gap in args.growth_gap:
                _apply_cost_scale(cost_scale)
                _apply_growth_gap(growth_gap)
                for death_floor in args.death_floor:
                    for selection_strength in selection_values:
                        config = replace(
                            base_config,
                            death_fitness_floor=death_floor,
                            selection_strength=selection_strength,
                        )

                        # Persistence test: no antibiotics.
                        mean_rf_off, _ = _resistant_fraction_trajectory(
                            config,
                            abx_on=False,
                            abx_class="none",
                            dose_level="std",
                            adherence=args.adherence,
                            immune_strength=args.immune_strength,
                            resistant_fraction=args.resistant_fraction,
                            dominant_genotype=args.dominant_genotype,
                            initial_population=args.initial_population,
                            n_days=args.n_days,
                            n_seeds=args.n_seeds,
                        )
                        half_life = _half_life_days(mean_rf_off)
                        rf_final = float(mean_rf_off[-1])
                        rf_day30 = float(mean_rf_off[min(29, len(mean_rf_off) - 1)])

                        # Counter-check: under beta-lactam, resistance should be
                        # selected (rf rises) and the population should not clear
                        # within 5 days.
                        mean_rf_on, clear_days_on = _resistant_fraction_trajectory(
                            config,
                            abx_on=True,
                            abx_class=args.counter_check_abx_class,
                            dose_level="std",
                            adherence=args.adherence,
                            immune_strength=args.immune_strength,
                            resistant_fraction=args.resistant_fraction,
                            dominant_genotype=args.dominant_genotype,
                            initial_population=args.initial_population,
                            n_days=min(30, args.n_days),
                            n_seeds=args.n_seeds,
                        )
                        cleared_early = [d for d in clear_days_on if d is not None and d <= 5]
                        rf_day5_on = float(mean_rf_on[min(4, len(mean_rf_on) - 1)])
                        rf_rises_under_abx = bool(rf_day5_on >= mean_rf_on[0])
                        counter_check_ok = rf_rises_under_abx and not cleared_early

                        meets_target = half_life >= args.target_half_life_days
                        grid_rows.append(
                            {
                                "cost_scale": cost_scale,
                                "growth_gap": growth_gap,
                                "death_fitness_floor": death_floor,
                                "selection_strength": selection_strength,
                                "rf_half_life_days": half_life,
                                "rf_final_no_abx": rf_final,
                                "rf_day30_no_abx": rf_day30,
                                "rf_day5_under_abx": rf_day5_on,
                                "n_cleared_early_under_abx": len(cleared_early),
                                "counter_check_ok": counter_check_ok,
                                "meets_persistence_target": meets_target,
                                "accept": bool(meets_target and counter_check_ok),
                            }
                        )
                        for day_index, value in enumerate(mean_rf_off, start=1):
                            traj_rows.append(
                                {
                                    "cost_scale": cost_scale,
                                    "growth_gap": growth_gap,
                                    "death_fitness_floor": death_floor,
                                    "selection_strength": selection_strength,
                                    "abx": "off",
                                    "day": day_index,
                                    "resistant_fraction": value,
                                }
                            )

                        flag = (
                            "ACCEPT"
                            if grid_rows[-1]["accept"]
                            else ("persist-ok" if meets_target else "washout")
                        )
                        print(
                            f"  cost={cost_scale:<4g} gap={growth_gap:<4g} "
                            f"floor={death_floor:<4g} sel={selection_strength:<4g} "
                            f"-> half_life={half_life:6.1f}d  rf_d30={rf_day30:.2f}  "
                            f"abx_d5={rf_day5_on:.2f}  "
                            f"check={'ok' if counter_check_ok else 'FAIL'}  [{flag}]"
                        )
    finally:
        _restore_levers()

    grid_df = pd.DataFrame(grid_rows)
    traj_df = pd.DataFrame(traj_rows)
    grid_path = data_dir / "micro_resistance_grid.parquet"
    traj_path = data_dir / "micro_resistance_trajectories.parquet"
    grid_df.to_parquet(grid_path, index=False)
    traj_df.to_parquet(traj_path, index=False)

    accepted = grid_df[grid_df["accept"]]
    print("=" * 72)
    if not accepted.empty:
        # Prefer the smallest deviation from current code among accepted points
        # (largest cost_scale, smallest death_floor, largest growth_gap), so the
        # recommended change stays minimal.
        best = accepted.sort_values(
            ["cost_scale", "growth_gap", "death_fitness_floor"],
            ascending=[False, False, True],
        ).iloc[0]
        print("Recommended setting (closest to current code that meets the target):")
        print(f"  ResistanceCosts scale = {best['cost_scale']:g}x  (genome.py)")
        print(
            f"  resistant GROWTH_BASE = {_SENSITIVE_GROWTH_BASE - best['growth_gap']:.2f}"
            f"  (engine.py::_create_seed_genome_for_genotype)"
        )
        print(f"  micro.death_fitness_floor = {best['death_fitness_floor']:g}  (YAML)")
        print(f"  micro.selection_strength  = {best['selection_strength']:g}  (YAML)")
        print(
            f"  -> rf_half_life = {best['rf_half_life_days']:g} days, "
            f"rf_day30 = {best['rf_day30_no_abx']:.2f}, "
            f"counter-check {'ok' if best['counter_check_ok'] else 'FAIL'}"
        )
    else:
        # Report the best-achievable so the user sees the ceiling explicitly.
        passed_check = grid_df[grid_df["counter_check_ok"]]
        pool = passed_check if not passed_check.empty else grid_df
        best = pool.sort_values("rf_half_life_days", ascending=False).iloc[0]
        print("No grid point met the persistence target AND the ABX counter-check.")
        print("This quantifies the abstraction ceiling: the residual fitness gap")
        print("between S and resistant strains is not fully reachable by these levers.")
        print("Best achievable (longest persistence with a passing ABX check):")
        print(
            f"  cost_scale={best['cost_scale']:g}  growth_gap={best['growth_gap']:g}  "
            f"death_fitness_floor={best['death_fitness_floor']:g}  "
            f"selection_strength={best['selection_strength']:g}"
        )
        print(
            f"  -> rf_half_life = {best['rf_half_life_days']:g} days, "
            f"rf_day30 = {best['rf_day30_no_abx']:.2f}, "
            f"counter-check {'ok' if best['counter_check_ok'] else 'FAIL'}"
        )
        print("Next lever if more persistence is needed: soften the resistant seed")
        print("genome's non-resistance genes (stress/dormancy/repair) in genome.py,")
        print("which carry implicit turnover costs beyond the three levers above.")
    print("\nSaved:")
    print(f"  Grid          -> {grid_path}")
    print(f"  Trajectories  -> {traj_path}")
    print("=" * 72)


if __name__ == "__main__":
    main()
