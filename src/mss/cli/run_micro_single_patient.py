"""Single-patient within-host micro simulation.

Runs one carrier episode through the full micro engine and generates
biology-focused plots of the within-host population dynamics.
No macro simulation — the only environmental input is the ABX schedule
defined in the config.

Outputs (outputs/<timestamp>_MicroSinglePatient/):
  data/single_patient_daily.parquet  -- daily metrics + mean gene values
  overview.png           -- total population, resistant fraction, diversity, p_clearance
  gene_expression.png    -- 14-gene heatmap: population-weighted mean trait over time
  genotype_composition.png -- stacked S/R1/R2/R3 fractions + total-pop overlay
  strain_landscape.png   -- block-max dominant strains as stacked fractions
  strain_landscape_total_population.png -- same lineages as absolute population
  lifecycle.png          -- mean damage load and mean lineage age

Usage:
    uv run python -m mss.cli.run_micro_single_patient
    uv run python -m mss.cli.run_micro_single_patient --config config/cal_micro_single_patient.yml
    uv run python -m mss.cli.run_micro_single_patient --output-dir outputs/my_run
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from mss.simulation.micro.single_patient import ABXPeriod, DayRecord, EpisodeConfig

matplotlib.use("Agg")

from mss.cli.run_coupled_simulation import PROJECT_ROOT
from mss.simulation.micro.config import SimulationConfig, build_micro_config
from mss.simulation.micro.engine import (
    StrainPopulation,
    compute_clearance_probability,
    simulate_day,
)
from mss.simulation.micro.genome import (
    NUM_GENES,
    classify_genotype,
    compute_resistant_fraction,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

GENE_NAMES: list[str] = [
    "Growth Base",
    "Metabolic Opt.",
    "Efflux Pumps",
    "Target Mod.",
    "Permeability",
    "Virulence",
    "Stealth",
    "Adhesion",
    "Mutation Rate",
    "HGT Competence",
    "DNA Repair",
    "Dormancy",
    "Stress Response",
    "Damage Tolerance",
]
assert len(GENE_NAMES) == NUM_GENES

# Gene group separators for the heatmap (index of first gene in next group)
_GENE_GROUP_BOUNDS = [2, 5, 8, 10]

GENO_COLORS: dict[str, str] = {
    "S": "#4c78a8",
    "R1": "#54a24b",
    "R2": "#f58518",
    "R3": "#e45756",
}
_ABX_SHADE = "#f28e2b"

# Per-genotype color palettes for the strain landscape (light→dark within class)
_STRAIN_PALETTES: dict[str, list[str]] = {
    "S": ["#4c78a8", "#6e9ec4", "#91bfe0", "#b4d9f5"],
    "R1": ["#54a24b", "#74bb6a", "#97d58c", "#bbecb0"],
    "R2": ["#f58518", "#f7a347", "#f9c176", "#fcdca5"],
    "R3": ["#e45756", "#ea7b7a", "#f0a0a0", "#f5c6c6"],
}

# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------


def _require_mapping(name: str, value: object) -> dict:
    if not isinstance(value, dict):
        raise TypeError(f"{name} must be a mapping.")
    return value


def _require_keys(name: str, raw: dict, keys: set[str]) -> None:
    missing = sorted(keys - set(raw))
    unknown = sorted(set(raw) - keys)
    if missing or unknown:
        parts: list[str] = []
        if missing:
            parts.append(f"missing keys: {', '.join(missing)}")
        if unknown:
            parts.append(f"unknown keys: {', '.join(unknown)}")
        raise ValueError(f"Invalid {name} config ({'; '.join(parts)}).")


def _load_config(path: Path) -> tuple[SimulationConfig, EpisodeConfig]:
    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict):
        raise TypeError("Top-level config must be a mapping.")

    micro_config = build_micro_config(_require_mapping("micro", raw.get("micro")), source="micro")

    ep_raw = _require_mapping("episode", raw.get("episode"))
    _require_keys(
        "episode",
        ep_raw,
        {
            "n_days",
            "seed",
            "resistant_fraction",
            "dominant_genotype",
            "initial_population",
            "immune_strength",
            "abx_schedule",
            "allow_spontaneous_clearance",
        },
    )
    schedule = []
    for index, item in enumerate(ep_raw["abx_schedule"]):
        period_raw = _require_mapping(f"episode.abx_schedule[{index}]", item)
        _require_keys(
            f"episode.abx_schedule[{index}]",
            period_raw,
            {"start_day", "end_day", "abx_class", "dose_level", "adherence"},
        )
        schedule.append(
            ABXPeriod(
                start_day=int(period_raw["start_day"]),
                end_day=int(period_raw["end_day"]),
                abx_class=str(period_raw["abx_class"]),
                dose_level=str(period_raw["dose_level"]),
                adherence=float(period_raw["adherence"]),
            )
        )

    episode_config = EpisodeConfig(
        n_days=int(ep_raw["n_days"]),
        seed=int(ep_raw["seed"]),
        resistant_fraction=float(ep_raw["resistant_fraction"]),
        dominant_genotype=str(ep_raw["dominant_genotype"]),
        initial_population=float(ep_raw["initial_population"]),
        immune_strength=float(ep_raw["immune_strength"]),
        abx_schedule=schedule,
        allow_spontaneous_clearance=bool(ep_raw.get("allow_spontaneous_clearance", True)),
    )
    return micro_config, episode_config


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------


def _abx_for_day(day: int, schedule: list[ABXPeriod]) -> tuple[str, str, float]:
    for period in schedule:
        if period.start_day <= day <= period.end_day:
            return period.abx_class, period.dose_level, period.adherence
    return "none", "std", 1.0


def _shannon_entropy(populations: np.ndarray) -> float:
    total = float(populations.sum())
    if total <= 0:
        return 0.0
    probs = populations / total
    probs = probs[probs > 0]
    return float(-np.sum(probs * np.log(probs)))


def run_episode(
    micro_config: SimulationConfig,
    ep_config: EpisodeConfig,
) -> list[DayRecord]:
    rng = np.random.default_rng(ep_config.seed)
    pop = StrainPopulation.create_initial(
        resistant_fraction=ep_config.resistant_fraction,
        dominant_genotype=ep_config.dominant_genotype,
        initial_population=ep_config.initial_population,
        rng=rng,
        strain_namespace="patient",
    )

    records: list[DayRecord] = []

    for day in range(1, ep_config.n_days + 1):
        abx_class, dose_level, adherence = _abx_for_day(day, ep_config.abx_schedule)

        pop, _ = simulate_day(
            population=pop,
            abx_class=abx_class,
            dose_level=dose_level,
            adherence=adherence,
            immune_strength=ep_config.immune_strength,
            config=micro_config,
            seed=ep_config.seed * 100_000 + day,
        )

        total_pop = float(pop.total_population)
        n = pop.genomes.shape[0]

        if total_pop > 0 and n > 0:
            weights = pop.populations / total_pop
            mean_genes = (weights[:, None] * pop.genomes).sum(axis=0).astype(np.float32)
            mean_damage = float((weights * pop.damage_loads).sum())
            mean_age = float((weights * pop.lineage_ages).sum())
        else:
            mean_genes = np.zeros(NUM_GENES, dtype=np.float32)
            mean_damage = 0.0
            mean_age = 0.0

        genotype_pop: dict[str, float] = {"S": 0.0, "R1": 0.0, "R2": 0.0, "R3": 0.0}
        strain_snapshot: list[tuple[str, float, str]] = []
        for i in range(n):
            geno = classify_genotype(pop.genomes[i])
            p = float(pop.populations[i])
            genotype_pop[geno] = genotype_pop.get(geno, 0.0) + p
            name = pop.strain_names[i] if pop.strain_names else f"strain_{i:04d}"
            strain_snapshot.append((name, p, geno))

        if total_pop > 0:
            frac = {k: v / total_pop for k, v in genotype_pop.items()}
        else:
            frac = {k: 0.0 for k in genotype_pop}

        p_clear = float(compute_clearance_probability(pop, ep_config.immune_strength, micro_config))

        records.append(
            DayRecord(
                day=day,
                total_population=total_pop,
                resistant_fraction=float(compute_resistant_fraction(pop.genomes, pop.populations)),
                n_strains=n,
                shannon_entropy=_shannon_entropy(pop.populations),
                p_clearance=p_clear,
                abx_class=abx_class,
                mean_genes=mean_genes,
                mean_damage=mean_damage,
                mean_age=mean_age,
                frac_S=frac.get("S", 0.0),
                frac_R1=frac.get("R1", 0.0),
                frac_R2=frac.get("R2", 0.0),
                frac_R3=frac.get("R3", 0.0),
                strain_snapshot=strain_snapshot,
                cleared=False,
            )
        )

        cleared_spontaneous = (
            rng.random() < p_clear if ep_config.allow_spontaneous_clearance else False
        )
        cleared = total_pop < micro_config.min_population or cleared_spontaneous

        if cleared:
            if cleared_spontaneous and total_pop >= micro_config.min_population:
                print(f"  Cleared spontaneously on day {day} (p_clearance={p_clear:.4f})")
            else:
                print(f"  Cleared on day {day}  (total_pop={total_pop:.1f})")
            records[-1].cleared = True
            break

    return records


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------


def _shade_abx(ax: plt.Axes, schedule: list[ABXPeriod], n_days: int) -> None:
    seen: set[str] = set()
    for p in schedule:
        label = f"ABX: {p.abx_class}" if p.abx_class not in seen else None
        ax.axvspan(
            max(0, p.start_day - 0.5),
            min(n_days + 0.5, p.end_day + 0.5),
            alpha=0.12,
            color=_ABX_SHADE,
            zorder=0,
            label=label,
        )
        seen.add(p.abx_class)


def _days(records: list[DayRecord]) -> np.ndarray:
    return np.array([r.day for r in records])


# ---------------------------------------------------------------------------
# Figure 1: Overview
# ---------------------------------------------------------------------------


def _fig_overview(records: list[DayRecord], ep: EpisodeConfig, out_dir: Path) -> None:
    days = _days(records)
    total_pop = np.array([r.total_population for r in records])
    rf = np.array([r.resistant_fraction for r in records])
    n_strains = np.array([r.n_strains for r in records])
    entropy = np.array([r.shannon_entropy for r in records])
    p_clear = np.array([r.p_clearance for r in records])

    fig, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
    fig.suptitle("Within-Host Dynamics — Overview", fontsize=14, fontweight="bold")

    # Total population (log scale)
    ax = axes[0, 0]
    ax.semilogy(days, np.clip(total_pop, 1, None), color="#333333", linewidth=1.8)
    ax.axhline(
        ep.initial_population,
        color="gray",
        linewidth=0.8,
        linestyle="--",
        alpha=0.5,
        label=f"initial  ({ep.initial_population:.0e})",
    )
    _shade_abx(ax, ep.abx_schedule, ep.n_days)
    ax.set_ylabel("Total bacterial population")
    ax.set_title("Total population")
    ax.legend(fontsize=8, loc="lower left")

    # Resistant fraction
    ax = axes[0, 1]
    ax.plot(days, rf, color=GENO_COLORS["R2"], linewidth=1.8)
    ax.axhline(0.5, color="gray", linewidth=0.8, linestyle="--", alpha=0.5)
    _shade_abx(ax, ep.abx_schedule, ep.n_days)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Resistant fraction")
    ax.set_title("Resistant fraction  (resistance score ≥ 0.3)")

    # Strain diversity + Shannon entropy
    ax = axes[1, 0]
    ax.plot(days, n_strains, color="#4c78a8", linewidth=1.8, label="Active strains")
    _shade_abx(ax, ep.abx_schedule, ep.n_days)
    ax2 = ax.twinx()
    ax2.plot(days, entropy, color="#9ecae1", linewidth=1.2, linestyle="--", label="Shannon entropy")
    ax2.set_ylabel("Shannon entropy", color="#9ecae1", fontsize=9)
    ax2.tick_params(axis="y", colors="#9ecae1")
    ax.set_ylabel("Active strains")
    ax.set_xlabel("Day")
    ax.set_title("Strain diversity")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper right")

    # p_clearance
    ax = axes[1, 1]
    ax.plot(days, p_clear, color="#54a24b", linewidth=1.8)
    ax.axhline(
        0.0039,
        color="gray",
        linewidth=0.8,
        linestyle="--",
        alpha=0.6,
        label="target ≈ 0.0039 / day",
    )
    _shade_abx(ax, ep.abx_schedule, ep.n_days)
    ax.set_ylim(0, 1)
    ax.set_ylabel("p_clearance / day")
    ax.set_xlabel("Day")
    ax.set_title("Daily clearance probability")
    ax.legend(fontsize=8)

    fig.tight_layout()
    path = out_dir / "overview.png"
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  {path.name}")


# ---------------------------------------------------------------------------
# Figure 2: Gene expression heatmap
# ---------------------------------------------------------------------------


def _fig_gene_expression(records: list[DayRecord], ep: EpisodeConfig, out_dir: Path) -> None:
    days = _days(records)
    # (NUM_GENES, n_days) — population-weighted mean trait value per gene per day
    gene_matrix = np.stack([r.mean_genes for r in records], axis=1)

    fig, ax = plt.subplots(figsize=(max(12, len(days) * 0.04 + 4), 6))
    im = ax.imshow(
        gene_matrix,
        aspect="auto",
        origin="upper",
        cmap="RdYlBu_r",
        vmin=0.0,
        vmax=1.0,
        extent=[days[0] - 0.5, days[-1] + 0.5, NUM_GENES - 0.5, -0.5],
    )
    ax.set_yticks(range(NUM_GENES))
    ax.set_yticklabels(GENE_NAMES, fontsize=9)
    ax.set_xlabel("Day")
    ax.set_title(
        "Mean gene traits across population  (population-weighted average)",
        fontweight="bold",
    )

    # Group separators (metabolism | resistance | virulence | evolvability | lifecycle)
    group_labels = ["Metabolism", "Resistance", "Virulence/Survival", "Evolvability", "Lifecycle"]
    group_starts = [0, 2, 5, 8, 10]
    for bound in _GENE_GROUP_BOUNDS:
        ax.axhline(bound - 0.5, color="white", linewidth=1.5, zorder=3)
    for label, start in zip(group_labels, group_starts):
        ax.text(
            days[0] - 0.5,
            start + 0.35,
            f"  {label}",
            fontsize=7,
            color="white",
            va="top",
            fontweight="bold",
            transform=ax.get_yaxis_transform(),
            ha="left",
        )

    # ABX periods as vertical bands
    for p in ep.abx_schedule:
        x0 = max(days[0] - 0.5, p.start_day - 0.5)
        x1 = min(days[-1] + 0.5, p.end_day + 0.5)
        ax.axvspan(x0, x1, alpha=0.22, color=_ABX_SHADE, zorder=2)
        ax.text(
            (x0 + x1) / 2,
            NUM_GENES - 0.4,
            p.abx_class.replace("_", "\n"),
            ha="center",
            va="bottom",
            fontsize=7,
            color="#a04000",
        )

    plt.colorbar(im, ax=ax, label="Mean trait value [0 – 1]", shrink=0.85, pad=0.01)
    fig.tight_layout()
    path = out_dir / "gene_expression.png"
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  {path.name}")


# ---------------------------------------------------------------------------
# Figure 3: Genotype composition (stacked area)
# ---------------------------------------------------------------------------


def _fig_genotype_composition(records: list[DayRecord], ep: EpisodeConfig, out_dir: Path) -> None:
    days = _days(records)
    total_pop = np.array([r.total_population for r in records])

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.stackplot(
        days,
        [r.frac_S for r in records],
        [r.frac_R1 for r in records],
        [r.frac_R2 for r in records],
        [r.frac_R3 for r in records],
        labels=["S – susceptible", "R1 – low resist.", "R2 – mid resist.", "R3 – high resist."],
        colors=[GENO_COLORS["S"], GENO_COLORS["R1"], GENO_COLORS["R2"], GENO_COLORS["R3"]],
        alpha=0.88,
    )
    _shade_abx(ax, ep.abx_schedule, ep.n_days)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Population fraction")
    ax.set_xlabel("Day")
    ax.set_title("Genotype composition over time", fontweight="bold")
    ax.legend(loc="upper left", fontsize=9)

    ax2 = ax.twinx()
    ax2.semilogy(
        days,
        np.clip(total_pop, 1, None),
        color="black",
        linewidth=1.0,
        linestyle=":",
        alpha=0.45,
        label="Total pop (log)",
    )
    ax2.set_ylabel("Total population (log)", color="gray", fontsize=9)
    ax2.tick_params(axis="y", colors="gray")
    ax2.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    path = out_dir / "genotype_composition.png"
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  {path.name}")


# ---------------------------------------------------------------------------
# Figure 4: Strain landscape
# ---------------------------------------------------------------------------


def _fig_strain_landscape(
    records: list[DayRecord],
    ep: EpisodeConfig,
    out_dir: Path,
    block_size_days: int = 100,
) -> None:
    days = _days(records)
    total_pop = np.array([r.total_population for r in records], dtype=np.float64)

    selected_strains: list[str] = []
    selected_set: set[str] = set()
    strain_peak_geno: dict[str, str] = {}
    block_edges: list[int] = []

    # Block maxima: select the single largest-fraction lineage in each time block.
    for block_start in range(int(days[0]), int(days[-1]) + 1, block_size_days):
        block_end = block_start + block_size_days
        block_edges.append(block_start)
        best_name = ""
        best_frac = -1.0
        best_geno = "S"
        for r in records:
            if not block_start <= r.day < block_end or r.total_population <= 0:
                continue
            for name, pop, geno in r.strain_snapshot:
                frac = pop / r.total_population
                if frac > best_frac:
                    best_name = name
                    best_frac = frac
                    best_geno = geno
        if best_name:
            strain_peak_geno[best_name] = best_geno
            if best_name not in selected_set:
                selected_set.add(best_name)
                selected_strains.append(best_name)

    name_idx = {name: i for i, name in enumerate(selected_strains)}
    fracs = np.zeros((len(selected_strains), len(days)), dtype=np.float64)
    abs_pops = np.zeros((len(selected_strains), len(days)), dtype=np.float64)
    for di, r in enumerate(records):
        if r.total_population <= 0:
            continue
        for name, pop, _ in r.strain_snapshot:
            if name in name_idx:
                idx = name_idx[name]
                abs_pops[idx, di] = pop
                fracs[idx, di] = pop / r.total_population

    other_frac = np.clip(1.0 - fracs.sum(axis=0), 0.0, 1.0)
    other_pop = np.clip(total_pop - abs_pops.sum(axis=0), 0.0, None)

    class_count: dict[str, int] = {}
    colors: list[str] = []
    for name in selected_strains:
        geno = strain_peak_geno.get(name, "S")
        palette = _STRAIN_PALETTES.get(geno, _STRAIN_PALETTES["S"])
        idx = class_count.get(geno, 0) % len(palette)
        colors.append(palette[idx])
        class_count[geno] = class_count.get(geno, 0) + 1

    legend_handles = [
        mpatches.Patch(facecolor=GENO_COLORS[g], label=f"{g} block-max lineages", alpha=0.88)
        for g in ["S", "R1", "R2", "R3"]
    ] + [mpatches.Patch(facecolor="#cccccc", label="other active strains", alpha=0.88)]

    def _draw_landscape(
        values: np.ndarray, other: np.ndarray, ylabel: str, title: str, filename: str
    ):
        fig, ax = plt.subplots(figsize=(14, 6))
        ax.stackplot(
            days,
            *values,
            other,
            colors=colors + ["#cccccc"],
            alpha=0.88,
        )
        for edge in block_edges[1:]:
            ax.axvline(edge - 0.5, color="white", linewidth=0.7, alpha=0.7, zorder=3)
        _shade_abx(ax, ep.abx_schedule, ep.n_days)
        ax.set_ylabel(ylabel)
        ax.set_xlabel("Day")
        ax.set_title(title, fontweight="bold")
        ax.legend(handles=legend_handles, loc="upper right", fontsize=9, framealpha=0.85)
        fig.tight_layout()
        path = out_dir / filename
        fig.savefig(path, dpi=130, bbox_inches="tight")
        plt.close(fig)
        print(f"  {path.name}")

    _draw_landscape(
        fracs,
        other_frac,
        "Population fraction",
        f"Strain landscape — block maxima, {block_size_days}-day blocks "
        f"({len(selected_strains)} unique lineages)",
        "strain_landscape.png",
    )
    _draw_landscape(
        abs_pops,
        other_pop,
        "Population size",
        f"Strain landscape — absolute population, {block_size_days}-day block maxima",
        "strain_landscape_total_population.png",
    )


# ---------------------------------------------------------------------------
# Figure 5: Lifecycle (damage + lineage age)
# ---------------------------------------------------------------------------


def _fig_lifecycle(records: list[DayRecord], ep: EpisodeConfig, out_dir: Path) -> None:
    days = _days(records)
    mean_damage = np.array([r.mean_damage for r in records])
    mean_age = np.array([r.mean_age for r in records])

    fig, axes = plt.subplots(2, 1, figsize=(14, 6), sharex=True)
    fig.suptitle("Lifecycle Dynamics", fontsize=13, fontweight="bold")

    ax = axes[0]
    ax.plot(days, mean_damage, color="#8b0000", linewidth=1.8)
    _shade_abx(ax, ep.abx_schedule, ep.n_days)
    ax.set_ylabel("Mean damage load")
    ax.set_title("Population-weighted mean accumulated damage")

    ax = axes[1]
    ax.plot(days, mean_age, color="#5b4ea8", linewidth=1.8)
    _shade_abx(ax, ep.abx_schedule, ep.n_days)
    ax.set_ylabel("Mean lineage age (steps)")
    ax.set_xlabel("Day")
    ax.set_title("Population-weighted mean lineage age")

    fig.tight_layout()
    path = out_dir / "lifecycle.png"
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  {path.name}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

_DEFAULT_CONFIG = PROJECT_ROOT / "config" / "cal_micro_single_patient.yml"


def main() -> None:
    parser = argparse.ArgumentParser(description="Single-patient within-host micro simulation.")
    parser.add_argument(
        "--config",
        type=Path,
        default=_DEFAULT_CONFIG,
        help=f"Config YAML (default: {_DEFAULT_CONFIG.name})",
    )
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    micro_config, ep_config = _load_config(args.config)

    output_dir = args.output_dir or (
        PROJECT_ROOT
        / "outputs"
        / (datetime.now().strftime("%Y%m%d_%H%M%S") + "_MicroSinglePatient")
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "data").mkdir(exist_ok=True)

    print("=" * 65)
    print("Single-patient within-host micro simulation")
    print(f"  Config:        {args.config.name}")
    print(f"  Duration:      {ep_config.n_days} days")
    print(
        f"  Initial state: {ep_config.resistant_fraction:.0%} resistant ({ep_config.dominant_genotype})"
    )
    print(f"  Immune str.:   {ep_config.immune_strength}")
    if ep_config.abx_schedule:
        for p in ep_config.abx_schedule:
            print(
                f"  ABX  days {p.start_day:>3}–{p.end_day:<3}  "
                f"{p.abx_class:<15}  {p.dose_level}  adherence={p.adherence}"
            )
    else:
        print("  ABX:  none")
    print("=" * 65)

    records = run_episode(micro_config, ep_config)

    cleared_day = records[-1].day if records[-1].cleared else None
    print(
        f"  Simulated {len(records)} days"
        + (f"  (cleared on day {cleared_day})" if cleared_day else "")
    )

    # Save daily data as parquet
    rows = []
    for r in records:
        row: dict = {
            "day": r.day,
            "total_population": r.total_population,
            "resistant_fraction": r.resistant_fraction,
            "n_strains": r.n_strains,
            "shannon_entropy": r.shannon_entropy,
            "p_clearance": r.p_clearance,
            "abx_class": r.abx_class,
            "mean_damage": r.mean_damage,
            "mean_age": r.mean_age,
            "frac_S": r.frac_S,
            "frac_R1": r.frac_R1,
            "frac_R2": r.frac_R2,
            "frac_R3": r.frac_R3,
        }
        for i, gname in enumerate(GENE_NAMES):
            key = f"gene_{gname.replace(' ', '_').replace('.', '').lower()}"
            row[key] = float(r.mean_genes[i])
        rows.append(row)

    pd.DataFrame(rows).to_parquet(output_dir / "data" / "single_patient_daily.parquet", index=False)

    print("\nGenerating figures...")
    _fig_overview(records, ep_config, output_dir)
    _fig_gene_expression(records, ep_config, output_dir)
    _fig_genotype_composition(records, ep_config, output_dir)
    _fig_strain_landscape(records, ep_config, output_dir, block_size_days=10)
    _fig_lifecycle(records, ep_config, output_dir)

    print("=" * 65)
    print(f"Output: {output_dir}")
    print("=" * 65)


if __name__ == "__main__":
    main()
