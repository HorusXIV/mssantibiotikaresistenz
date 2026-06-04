"""Visualize simulation results from a run's data directory.

Usage:
    python visualize_results.py
    python visualize_results.py --data-dir outputs/runs/<timestamp>/data --plot-dir outputs/runs/<timestamp>/plots
    python visualize_results.py --interactive

Produces diagnostic plots that show whether the simulation is behaving
realistically:
  01_population_overview.png   — population size, prevalence, resistance, clinical load
  02_occupancy_stability.png   — per-hospital occupancy (validates Poisson ↔ discharge)
  03_hospital_prevalence.png   — per-hospital carrier prevalence (spatial distribution)
  04_department_mix.png        — Ward / ICU patient counts + isolated patients over time
  05_hospital_grid_heatmap.png — 3×2 grid coloured by final-day carrier prevalence
  06_department_grid.png       — per-hospital department-zone density view
  07_micro_evolution_overview.png  — micro dynamics over time
  08_micro_genotype_stream.png     — dominant-genotype composition over time
  09_micro_hospital_heatmap.png    — hospital-level mean resistance over time
  10_micro_patient_phase_space.png — patient-level resistance vs clearance (final day)
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.widgets import Slider

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_DATA_DIR = PROJECT_ROOT / "outputs" / "data"
DEFAULT_PLOT_DIR = PROJECT_ROOT / "outputs" / "plots"

# Expected endemic prevalence band (for reference lines)
ENDEMIC_LOW = 0.15
ENDEMIC_HIGH = 0.25

ROLLING_WINDOW = 7  # days for smoothing


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _load(
    data_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    daily_path = data_dir / "macro_daily.parquet"
    hosp_path = data_dir / "macro_daily_by_hospital.parquet"
    micro_daily_path = data_dir / "micro_daily.parquet"
    micro_hosp_path = data_dir / "micro_daily_by_hospital.parquet"
    micro_patient_path = data_dir / "micro_patient_daily.parquet"
    micro_genotype_path = data_dir / "micro_daily_genotype.parquet"
    if not daily_path.exists():
        raise FileNotFoundError(f"Parquet not found: {daily_path}\nRun the simulation first.")
    daily = pd.read_parquet(daily_path)
    hosp = pd.read_parquet(hosp_path) if hosp_path.exists() else pd.DataFrame()
    micro_daily = pd.read_parquet(micro_daily_path) if micro_daily_path.exists() else pd.DataFrame()
    micro_hosp = pd.read_parquet(micro_hosp_path) if micro_hosp_path.exists() else pd.DataFrame()
    micro_patient = (
        pd.read_parquet(micro_patient_path) if micro_patient_path.exists() else pd.DataFrame()
    )
    micro_genotype = (
        pd.read_parquet(micro_genotype_path) if micro_genotype_path.exists() else pd.DataFrame()
    )
    return daily, hosp, micro_daily, micro_hosp, micro_patient, micro_genotype


_quiet_mode = False  # set by run() before plotting


def _save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    if not _quiet_mode:
        print(f"  saved → {path}")


def _rolling(series: pd.Series, window: int = ROLLING_WINDOW) -> pd.Series:
    return series.rolling(window, min_periods=1, center=True).mean()


def _series(df: pd.DataFrame, column: str, default: float = 0.0) -> pd.Series:
    if column in df.columns:
        return pd.to_numeric(df[column], errors="coerce").fillna(default)
    return pd.Series(default, index=df.index, dtype=float)


# ---------------------------------------------------------------------------
# Plot 1: Population Overview (2×2)
# ---------------------------------------------------------------------------


def plot_population_overview(daily: pd.DataFrame, plot_dir: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    fig.suptitle("Population Overview", fontsize=14, fontweight="bold")
    days = daily["day"]

    # --- [0,0] Absolute counts ---
    ax = axes[0, 0]
    ax.plot(days, daily["total_patients"], color="gray", lw=1.5, label="Total")
    ax.plot(days, daily["susceptible"], color="steelblue", lw=1.5, label="Susceptible")
    ax.plot(days, daily["carriers"], color="firebrick", lw=1.5, label="Carriers")
    ax.set_title("Patient Counts")
    ax.set_xlabel("Day (d)")
    ax.set_ylabel("Patients (count)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # --- [0,1] Carrier prevalence ---
    ax = axes[0, 1]
    prev = daily["prevalence"] * 100
    ax.plot(days, prev, color="firebrick", lw=1, alpha=0.5, label="Daily")
    ax.plot(days, _rolling(prev), color="firebrick", lw=2, label=f"{ROLLING_WINDOW}d mean")
    ax.set_title("Carrier Prevalence (%)")
    ax.set_xlabel("Day (d)")
    ax.set_ylabel("Prevalence (%)")
    ax.set_ylim(0, max(prev.max() * 1.15, 30))
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # --- [1,0] Resistance evolution ---
    ax = axes[1, 0]
    res = daily["avg_resistant_fraction"]
    ax.plot(days, res, color="darkorange", lw=1, alpha=0.5, label="Daily")
    ax.plot(days, _rolling(res), color="darkorange", lw=2, label=f"{ROLLING_WINDOW}d mean")
    ax.set_title("Avg Resistant Fraction (carriers, share)")
    ax.set_xlabel("Day (d)")
    ax.set_ylabel("Resistant fraction (share)")
    ax.set_ylim(0, 1)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # --- [1,1] Clinical load ---
    ax = axes[1, 1]
    ax.plot(days, daily["isolated_count"], color="purple", lw=1.5, label="Isolated (flag)")
    ax.plot(days, daily["abx_on_count"], color="teal", lw=1.5, label="On ABX")
    ax.set_title("Clinical Load (Isolation flag & ABX)")
    ax.set_xlabel("Day (d)")
    ax.set_ylabel("Patients (count)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    _save(fig, plot_dir / "01_population_overview.png")


# ---------------------------------------------------------------------------
# Plot 2: Occupancy Stability
# ---------------------------------------------------------------------------


def plot_occupancy_stability(hosp: pd.DataFrame, plot_dir: Path) -> None:
    if hosp.empty or "hospital_id" not in hosp.columns:
        print("  skipped 02_occupancy_stability.png (no per-hospital data)")
        return

    fig, ax = plt.subplots(figsize=(11, 5))
    fig.suptitle("Hospital Occupancy (patients)", fontsize=14, fontweight="bold")

    hospital_ids = sorted(hosp["hospital_id"].unique())
    colors = plt.cm.tab10(np.linspace(0, 0.9, len(hospital_ids)))

    for hid, color in zip(hospital_ids, colors):
        sub = hosp[hosp["hospital_id"] == hid].sort_values("day")
        ax.plot(
            sub["day"],
            sub["total_patients"],
            color=color,
            lw=1.2,
            alpha=0.8,
            label=hid.replace("hospital_", "H"),
        )

    # Global mean
    global_mean = hosp.groupby("day")["total_patients"].mean()
    ax.plot(
        global_mean.index,
        global_mean.values,
        color="black",
        lw=2,
        linestyle="--",
        label="Global mean",
        zorder=5,
    )

    ax.set_xlabel("Day (d)")
    ax.set_ylabel("Patients (count)")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3)

    # Annotate final values
    last_day = hosp["day"].max()
    for hid, color in zip(hospital_ids, colors):
        sub = hosp[(hosp["hospital_id"] == hid) & (hosp["day"] == last_day)]
        if not sub.empty:
            val = sub["total_patients"].values[0]
            ax.annotate(f"{val:.0f}", xy=(last_day, val), fontsize=7, color=color, va="center")

    fig.tight_layout()
    _save(fig, plot_dir / "02_occupancy_stability.png")


# ---------------------------------------------------------------------------
# Plot 3: Per-Hospital Carrier Prevalence
# ---------------------------------------------------------------------------


def plot_hospital_prevalence(hosp: pd.DataFrame, plot_dir: Path) -> None:
    if hosp.empty:
        print("  skipped 03_hospital_prevalence.png (no per-hospital data)")
        return

    fig, ax = plt.subplots(figsize=(11, 5))
    fig.suptitle("Carrier Prevalence by Hospital (%)", fontsize=14, fontweight="bold")

    hospital_ids = sorted(hosp["hospital_id"].unique())
    colors = plt.cm.tab10(np.linspace(0, 0.9, len(hospital_ids)))

    for hid, color in zip(hospital_ids, colors):
        sub = hosp[hosp["hospital_id"] == hid].sort_values("day")
        prev = sub["prevalence"] * 100
        ax.plot(sub["day"], prev, color=color, lw=1.2, label=hid.replace("hospital_", "H"))

    # Global mean
    global_prev = hosp.groupby("day").apply(
        lambda g: (
            g["carriers"].sum() / g["total_patients"].sum() * 100
            if g["total_patients"].sum() > 0
            else 0
        ),
        include_groups=False,
    )
    ax.plot(
        global_prev.index,
        global_prev.values,
        color="black",
        lw=2.5,
        linestyle="--",
        label="Global mean",
        zorder=5,
    )

    ax.set_xlabel("Day (d)")
    ax.set_ylabel("Carrier prevalence (%)")
    ax.set_ylim(0, None)
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    _save(fig, plot_dir / "03_hospital_prevalence.png")


# ---------------------------------------------------------------------------
# Plot 4: Department Mix
# ---------------------------------------------------------------------------


def plot_department_mix(daily: pd.DataFrame, plot_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 5))
    fig.suptitle("Department Mix (patients)", fontsize=14, fontweight="bold")

    days = daily["day"].values
    ward = daily["ward_count"].values
    icu = daily["icu_count"].values

    ax.stackplot(
        days, ward, icu, labels=["Ward", "ICU"], colors=["steelblue", "firebrick"], alpha=0.65
    )

    # Isolated patients: is_isolated flag (not a physical grid zone — carriers
    # detected and flagged; reduces their transmission by iso_reduction factor)
    if "isolated_count" in daily.columns:
        ax.plot(
            days,
            daily["isolated_count"].values,
            color="purple",
            lw=2,
            label="Isolated (detection flag)",
            zorder=5,
        )

    ax.set_xlabel("Day (d)")
    ax.set_ylabel("Patients (count)")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(alpha=0.2)

    fig.tight_layout()
    _save(fig, plot_dir / "04_department_mix.png")


# ---------------------------------------------------------------------------
# Plot 5: Hospital Grid Heatmap (final day)
# ---------------------------------------------------------------------------


def plot_hospital_grid_heatmap(hosp: pd.DataFrame, plot_dir: Path, grid_cols: int = 3) -> None:
    if hosp.empty:
        print("  skipped 05_hospital_grid_heatmap.png (no per-hospital data)")
        return

    last_day = hosp["day"].max()
    final = hosp[hosp["day"] == last_day].copy()
    final = final[final["total_patients"] > 0].copy()
    final["prev_pct"] = final["carriers"] / final["total_patients"] * 100

    hospital_ids = sorted(final["hospital_id"].unique())
    n = len(hospital_ids)
    grid_rows = math.ceil(n / grid_cols)

    # Build grid matrix
    prev_grid = np.full((grid_rows, grid_cols), np.nan)
    occ_grid = np.zeros((grid_rows, grid_cols))
    label_grid: list[list[str]] = [[""] * grid_cols for _ in range(grid_rows)]

    for hid in hospital_ids:
        idx = int(hid.split("_")[-1]) - 1
        row, col = divmod(idx, grid_cols)
        if row < grid_rows:
            row_data = final[final["hospital_id"] == hid]
            if not row_data.empty:
                prev_grid[row, col] = row_data["prev_pct"].values[0]
                occ_grid[row, col] = row_data["total_patients"].values[0]
                label_grid[row][col] = (
                    f"{hid.replace('hospital_', 'H')}\n"
                    f"prev: {prev_grid[row, col]:.1f}%\n"
                    f"n={occ_grid[row, col]:.0f}"
                )

    fig, ax = plt.subplots(figsize=(8, 4))
    fig.suptitle(
        f"Hospital Grid — Carrier Prevalence (%) (day {last_day})",
        fontsize=13,
        fontweight="bold",
    )

    im = ax.imshow(
        prev_grid,
        cmap="RdYlGn_r",
        vmin=0,
        vmax=max(ENDEMIC_HIGH * 100 * 1.5, np.nanmax(prev_grid) * 1.05),
    )

    for r in range(grid_rows):
        for c in range(grid_cols):
            if not np.isnan(prev_grid[r, c]):
                text_color = "white" if prev_grid[r, c] > ENDEMIC_HIGH * 100 * 1.2 else "black"
                ax.text(
                    c, r, label_grid[r][c], ha="center", va="center", fontsize=9, color=text_color
                )

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label("Carrier prevalence (%)")

    ax.set_xticks(range(grid_cols))
    ax.set_yticks(range(grid_rows))
    ax.set_xticklabels([f"Col {c}" for c in range(grid_cols)])
    ax.set_yticklabels([f"Row {r}" for r in range(grid_rows)])
    ax.set_title("Low -> High prevalence (%)", fontsize=9)

    fig.tight_layout()
    _save(fig, plot_dir / "05_hospital_grid_heatmap.png")


# ---------------------------------------------------------------------------
# Plot 6: Department Grid per Hospital (final day)
# ---------------------------------------------------------------------------


def plot_department_grid(cell_df: pd.DataFrame, plot_dir: Path) -> None:
    """Per-cell carrier prevalence heatmap for every hospital (final day).

    Each grid cell is coloured by its carrier prevalence (carriers / patients).
    This is the spatial structure proximity_decay_alpha actually acts on: a
    susceptible's infection risk is dominated by carriers in its own and
    neighbouring cells. Grid dimensions and the ICU/Ward split are read from the
    logged data, so the plot always matches the simulated grid. Each cell is
    annotated with its prevalence and patient count n.
    """
    if cell_df.empty:
        print("  skipped 06_department_grid.png (no per-cell data)")
        return

    last_day = cell_df["day"].max()
    final = cell_df[cell_df["day"] == last_day].copy()
    hospital_ids = sorted(final["hospital_id"].unique())
    n_hospitals = len(hospital_ids)

    dept_cols = int(final["x"].max()) + 1
    dept_rows = int(final["y"].max()) + 1
    vmax = max(5.0, float(final["prevalence"].max()) * 100.0)

    ncols_fig = min(n_hospitals, 3)
    nrows_fig = math.ceil(n_hospitals / ncols_fig)
    fig, axes = plt.subplots(
        nrows_fig,
        ncols_fig,
        figsize=(4.5 * ncols_fig, 3.2 * nrows_fig),
        constrained_layout=True,
    )
    fig.suptitle(
        "Department Grid — Carrier Prevalence per Cell (final day)",
        fontsize=13,
        fontweight="bold",
    )
    axes = np.array(axes).reshape(nrows_fig, ncols_fig)

    im = None
    for idx, hid in enumerate(hospital_ids):
        r, c = divmod(idx, ncols_fig)
        ax = axes[r, c]
        sub = final[final["hospital_id"] == hid]
        if sub.empty:
            ax.axis("off")
            continue

        prev_grid = np.zeros((dept_rows, dept_cols))
        labels: list[list[str]] = [[""] * dept_cols for _ in range(dept_rows)]
        row_dept: dict[int, str] = {}
        for _, cell in sub.iterrows():
            x, y = int(cell["x"]), int(cell["y"])
            prev_grid[y, x] = float(cell["prevalence"]) * 100.0
            labels[y][
                x
            ] = f"{float(cell['prevalence']) * 100:.0f}%\nn={int(cell['total_patients'])}"
            row_dept[y] = str(cell["department"])

        im = ax.imshow(prev_grid, cmap="YlOrRd", vmin=0.0, vmax=vmax, aspect="auto")
        for y in range(dept_rows):
            for x in range(dept_cols):
                ax.text(x, y, labels[y][x], ha="center", va="center", fontsize=7, color="black")

        # Separate ICU rows (top) from Ward rows with a line.
        icu_ys = [y for y, d in row_dept.items() if d == "icu"]
        if icu_ys and len(icu_ys) < dept_rows:
            ax.axhline(max(icu_ys) + 0.5, color="black", lw=1.5)

        ax.set_title(hid.replace("hospital_", "H"), fontsize=10)
        ax.set_xticks(range(dept_cols))
        ax.set_yticks(range(dept_rows))
        ax.set_xticklabels([f"C{i}" for i in range(dept_cols)], fontsize=7)
        ax.set_yticklabels([row_dept.get(y, "")[:4] for y in range(dept_rows)], fontsize=7)

    for idx in range(n_hospitals, nrows_fig * ncols_fig):
        r, c = divmod(idx, ncols_fig)
        axes[r, c].axis("off")

    if im is not None:
        fig.colorbar(
            im,
            ax=axes.ravel().tolist(),
            fraction=0.025,
            pad=0.02,
            label="Carrier prevalence (%)",
        )
    _save(fig, plot_dir / "06_department_grid.png")


# ---------------------------------------------------------------------------
# Plot 7: Micro Evolution Overview
# ---------------------------------------------------------------------------


def plot_micro_evolution_overview(micro_daily: pd.DataFrame, plot_dir: Path) -> None:
    if micro_daily.empty or "day" not in micro_daily.columns:
        print("  skipped 07_micro_evolution_overview.png (no micro daily data)")
        return

    micro_daily = micro_daily.sort_values("day")
    days = micro_daily["day"]

    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    fig.suptitle("Micro Evolution Overview", fontsize=14, fontweight="bold")

    carriers = _series(micro_daily, "carrier_count")
    active = _series(micro_daily, "active_episodes")
    ax = axes[0, 0]
    ax.plot(days, carriers, lw=1.8, color="steelblue", label="Carrier count")
    ax.plot(days, active, lw=1.8, color="darkorange", label="Active micro episodes")
    ax.set_title("Carrier Load")
    ax.set_xlabel("Day (d)")
    ax.set_ylabel("Count (patients/episodes)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    p10 = _series(micro_daily, "p10_resistant_fraction") * 100
    p50 = _series(micro_daily, "p50_resistant_fraction") * 100
    p90 = _series(micro_daily, "p90_resistant_fraction") * 100
    ax = axes[0, 1]
    ax.fill_between(days, p10, p90, color="firebrick", alpha=0.2, label="P10-P90")
    ax.plot(days, p50, color="firebrick", lw=2, label="Median")
    ax.plot(days, _rolling(p50), color="black", lw=1.2, ls="--", label=f"{ROLLING_WINDOW}d mean")
    ax.set_title("Within-Host Resistance Distribution (%)")
    ax.set_xlabel("Day (d)")
    ax.set_ylabel("Resistant fraction (%)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    mean_n_strains = _series(micro_daily, "mean_n_strains")
    mean_total_pop = _series(micro_daily, "mean_total_population")
    ax = axes[1, 0]
    ln1 = ax.plot(days, mean_n_strains, color="purple", lw=2, label="Mean strains / episode")
    ax.set_xlabel("Day (d)")
    ax.set_ylabel("Strains (count)", color="purple")
    ax.tick_params(axis="y", labelcolor="purple")
    ax.grid(alpha=0.3)
    ax2 = ax.twinx()
    ln2 = ax2.plot(
        days,
        np.log10(mean_total_pop + 1.0),
        color="teal",
        lw=1.8,
        label=r"log10(mean total population + 1)",
    )
    ax2.set_ylabel(r"log10 population (a.u.)", color="teal")
    ax2.tick_params(axis="y", labelcolor="teal")
    ax.set_title("Strain Diversity and Bacterial Load")
    lines = ln1 + ln2
    ax.legend(lines, [line.get_label() for line in lines], fontsize=8, loc="upper left")

    p_clear = _series(micro_daily, "mean_p_clearance") * 100
    transmissibility = _series(micro_daily, "mean_relative_transmissibility")
    ax = axes[1, 1]
    ax.plot(days, p_clear, color="green", lw=2, label="Mean p_clearance (%)")
    ax.plot(
        days,
        transmissibility,
        color="darkred",
        lw=2,
        label="Mean relative transmissibility",
    )
    ax.set_title("Micro-Derived Clinical Effectors")
    ax.set_xlabel("Day (d)")
    ax.set_ylabel("Value (%, rel.)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    _save(fig, plot_dir / "07_micro_evolution_overview.png")


# ---------------------------------------------------------------------------
# Plot 8: Micro Dominant-Genotype Stream
# ---------------------------------------------------------------------------


def plot_micro_genotype_stream(micro_genotype: pd.DataFrame, plot_dir: Path) -> None:
    if micro_genotype.empty or "day" not in micro_genotype.columns:
        print("  skipped 08_micro_genotype_stream.png (no micro genotype data)")
        return

    pivot = (
        micro_genotype.pivot_table(
            index="day",
            columns="dominant_genotype",
            values="fraction",
            aggfunc="sum",
            fill_value=0.0,
        )
        .sort_index()
        .fillna(0.0)
    )
    if pivot.empty:
        print("  skipped 08_micro_genotype_stream.png (empty genotype pivot)")
        return

    preferred_order = ["S", "R1", "R2", "R3", "OTHER"]
    columns = [c for c in preferred_order if c in pivot.columns] + [
        c for c in pivot.columns if c not in preferred_order
    ]
    pivot = pivot[columns]

    color_map = {
        "S": "#4daf4a",
        "R1": "#ff7f00",
        "R2": "#e41a1c",
        "R3": "#984ea3",
        "OTHER": "#7f7f7f",
    }
    colors = [color_map.get(col, "#7f7f7f") for col in pivot.columns]

    fig, ax = plt.subplots(figsize=(11, 5))
    fig.suptitle("Dominant Genotype Composition Over Time", fontsize=14, fontweight="bold")
    ax.stackplot(
        pivot.index, [pivot[c] * 100 for c in pivot.columns], labels=pivot.columns, colors=colors
    )
    ax.set_xlabel("Day (d)")
    ax.set_ylabel("Share of carriers (%)")
    ax.set_ylim(0, 100)
    ax.legend(loc="upper left", fontsize=8, ncol=3)
    ax.grid(alpha=0.25)

    fig.tight_layout()
    _save(fig, plot_dir / "08_micro_genotype_stream.png")


# ---------------------------------------------------------------------------
# Plot 9: Micro Hospital Resistance Heatmap (day × hospital)
# ---------------------------------------------------------------------------


def plot_micro_hospital_heatmap(micro_hosp: pd.DataFrame, plot_dir: Path) -> None:
    if micro_hosp.empty or "hospital_id" not in micro_hosp.columns:
        print("  skipped 09_micro_hospital_heatmap.png (no micro hospital data)")
        return

    heat = (
        micro_hosp.assign(
            mean_resistant_fraction_pct=_series(micro_hosp, "mean_resistant_fraction") * 100
        )
        .pivot_table(
            index="hospital_id",
            columns="day",
            values="mean_resistant_fraction_pct",
            aggfunc="mean",
        )
        .sort_index()
        .fillna(0.0)
    )
    if heat.empty:
        print("  skipped 09_micro_hospital_heatmap.png (empty heatmap data)")
        return

    fig, ax = plt.subplots(figsize=(12, max(4, 0.5 * len(heat.index) + 2)))
    fig.suptitle("Hospital Resistance Heatmap (%)", fontsize=14, fontweight="bold")
    im = ax.imshow(
        heat.values, aspect="auto", cmap="YlOrRd", vmin=0, vmax=max(25.0, np.max(heat.values))
    )
    ax.set_xlabel("Day (d)")
    ax.set_ylabel("Hospital")
    ax.set_xticks(range(len(heat.columns)))
    ax.set_xticklabels([str(int(day)) for day in heat.columns], rotation=90, fontsize=7)
    ax.set_yticks(range(len(heat.index)))
    ax.set_yticklabels([h.replace("hospital_", "H") for h in heat.index], fontsize=8)
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Mean resistant fraction (%)")

    fig.tight_layout()
    _save(fig, plot_dir / "09_micro_hospital_heatmap.png")


# ---------------------------------------------------------------------------
# Plot 10: Micro Patient Phase Space (final day)
# ---------------------------------------------------------------------------


def plot_micro_patient_phase_space(micro_patient: pd.DataFrame, plot_dir: Path) -> None:
    if micro_patient.empty or "day" not in micro_patient.columns:
        print("  skipped 10_micro_patient_phase_space.png (no micro patient data)")
        return

    last_day = int(micro_patient["day"].max())
    final = micro_patient[micro_patient["day"] == last_day].copy()
    if final.empty:
        print("  skipped 10_micro_patient_phase_space.png (empty final-day micro data)")
        return

    x = _series(final, "resistant_fraction") * 100
    y = _series(final, "p_clearance") * 100
    size = 18 + 6 * np.clip(_series(final, "n_strains"), 0, 25)
    genotype = final["dominant_genotype"].astype(str)
    color_map = {"S": "#4daf4a", "R1": "#ff7f00", "R2": "#e41a1c", "R3": "#984ea3"}
    colors = [color_map.get(g, "#7f7f7f") for g in genotype]

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.suptitle(f"Micro Patient Phase Space (day {last_day})", fontsize=14, fontweight="bold")
    ax.scatter(x, y, s=size, c=colors, alpha=0.65, edgecolor="black", linewidth=0.2)
    ax.axvline(np.median(x) if len(x) > 0 else 0, ls="--", lw=1, color="gray", alpha=0.5)
    ax.axhline(np.median(y) if len(y) > 0 else 0, ls="--", lw=1, color="gray", alpha=0.5)
    ax.set_xlabel("Resistant fraction (%)")
    ax.set_ylabel("p_clearance (%)")
    ax.grid(alpha=0.3)

    legend_labels = ["S", "R1", "R2", "R3", "OTHER"]
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=color_map.get(label, "#7f7f7f"),
            markeredgecolor="black",
            markeredgewidth=0.2,
            markersize=7,
            label=label,
        )
        for label in legend_labels
    ]
    ax.legend(handles=handles, title="Genotype", fontsize=8, title_fontsize=9, loc="upper right")

    fig.tight_layout()
    _save(fig, plot_dir / "10_micro_patient_phase_space.png")


# ---------------------------------------------------------------------------
# Optional interactive explorer (micro patient points by day)
# ---------------------------------------------------------------------------


def launch_micro_patient_explorer(micro_patient: pd.DataFrame) -> None:
    if micro_patient.empty or "day" not in micro_patient.columns:
        print("  skipped interactive micro explorer (no micro patient data)")
        return

    days = sorted(int(d) for d in micro_patient["day"].dropna().unique())
    if not days:
        print("  skipped interactive micro explorer (no valid day values)")
        return

    fig, ax = plt.subplots(figsize=(11, 6))
    fig.subplots_adjust(bottom=0.18)

    slider_ax = fig.add_axes([0.12, 0.06, 0.76, 0.04])
    slider = Slider(
        ax=slider_ax,
        label="Day",
        valmin=min(days),
        valmax=max(days),
        valinit=days[-1],
        valstep=days,
    )

    color_map = {"S": "#4daf4a", "R1": "#ff7f00", "R2": "#e41a1c", "R3": "#984ea3"}
    scatter = ax.scatter([], [], s=[], alpha=0.65, edgecolor="black", linewidth=0.2)

    def _update(day_value: float) -> None:
        day = int(day_value)
        sub = micro_patient[micro_patient["day"] == day]
        x = _series(sub, "resistant_fraction") * 100
        y = _series(sub, "p_clearance") * 100
        sizes = 18 + 6 * np.clip(_series(sub, "n_strains"), 0, 25)
        colors = [color_map.get(str(g), "#7f7f7f") for g in sub["dominant_genotype"].astype(str)]

        points = np.column_stack([x.values, y.values]) if len(sub) > 0 else np.empty((0, 2))
        scatter.set_offsets(points)
        scatter.set_sizes(sizes.values if len(sub) > 0 else [])
        scatter.set_color(colors)

        ax.set_xlim(0, max(100.0, float(x.max()) * 1.1 if len(sub) > 0 else 100.0))
        ax.set_ylim(0, max(100.0, float(y.max()) * 1.2 if len(sub) > 0 else 100.0))
        ax.set_xlabel("Resistant fraction (%)")
        ax.set_ylabel("p_clearance (%)")
        ax.set_title(f"Interactive Micro Patient Explorer — day {day} (n={len(sub)})")
        ax.grid(alpha=0.3)
        fig.canvas.draw_idle()

    slider.on_changed(_update)
    _update(days[-1])
    plt.show()


# ---------------------------------------------------------------------------
# Sanity-check summary printed to console
# ---------------------------------------------------------------------------


def _print_sanity_checks(
    daily: pd.DataFrame,
    hosp: pd.DataFrame,
    micro_daily: pd.DataFrame,
) -> None:
    last = daily.iloc[-1]
    final_prev = last["prevalence"] * 100
    final_res = last["avg_resistant_fraction"]
    min_occ = daily["total_patients"].min()
    max_occ = daily["total_patients"].max()
    mean_occ = daily["total_patients"].mean()

    print("\n── Sanity Checks ──────────────────────────────────────────")
    print(f"  Days simulated:        {int(last['day'])}")
    print(f"  Final prevalence:      {final_prev:.1f}%")
    print(f"  Final avg resistance:  {final_res:.3f}")
    print(f"  Occupancy (min/mean/max): {min_occ:.0f} / {mean_occ:.0f} / {max_occ:.0f}")

    if not hosp.empty:
        last_hosp = hosp[hosp["day"] == hosp["day"].max()]
        per_h = last_hosp.set_index("hospital_id")["total_patients"]
        print(f"  Final occupancy range: {per_h.min():.0f}–{per_h.max():.0f} per hospital")
        spread = per_h.max() - per_h.min()
        print(
            f"  Inter-hospital spread: {spread:.0f} patients"
            + (" (balanced)" if spread < per_h.mean() * 0.5 else " (unbalanced)")
        )

    occ_cv = daily["total_patients"].std() / daily["total_patients"].mean()
    print(
        f"  Occupancy CV:          {occ_cv:.3f}"
        + (" (stable)" if occ_cv < 0.10 else " (volatile — check admission rate)")
    )
    if not micro_daily.empty:
        m_last = micro_daily.iloc[-1]
        print(f"  Final active episodes: {int(m_last.get('active_episodes', 0))}")
        print(f"  Final mean p_clear:    {float(m_last.get('mean_p_clearance', 0.0)):.3f}")
        print(f"  Final mean n_strains:  {float(m_last.get('mean_n_strains', 0.0)):.2f}")
    print("────────────────────────────────────────────────────────────\n")


# ---------------------------------------------------------------------------
# Core entry point (importable)
# ---------------------------------------------------------------------------


def run(
    data_dir: Path = DEFAULT_DATA_DIR,
    plot_dir: Path = DEFAULT_PLOT_DIR,
    quiet: bool = False,
    interactive: bool = False,
) -> None:
    """Generate all diagnostic plots from simulation parquet files.

    Parameters
    ----------
    quiet
        When True (used by run_coupled_simulation.py), suppress all output
        except a single summary line printed at the end.
    """
    global _quiet_mode
    _quiet_mode = quiet

    daily, hosp, micro_daily, micro_hosp, micro_patient, micro_genotype = _load(data_dir)

    cell_path = data_dir / "macro_cell_daily.parquet"
    cell_df = pd.read_parquet(cell_path) if cell_path.exists() else pd.DataFrame()

    if not quiet:
        print(f"Reading data from: {data_dir}")
        print(f"  macro_daily:             {len(daily)} rows")
        print(f"  macro_daily_by_hospital: {len(hosp)} rows")
        print(f"  micro_daily:             {len(micro_daily)} rows")
        print(f"  micro_daily_by_hospital: {len(micro_hosp)} rows")
        print(f"  micro_patient_daily:     {len(micro_patient)} rows")
        print(f"  micro_daily_genotype:    {len(micro_genotype)} rows")
        print(f"\nGenerating plots → {plot_dir}")

    plot_population_overview(daily, plot_dir)
    plot_occupancy_stability(hosp, plot_dir)
    plot_hospital_prevalence(hosp, plot_dir)
    plot_department_mix(daily, plot_dir)
    plot_hospital_grid_heatmap(hosp, plot_dir)
    plot_department_grid(cell_df, plot_dir)
    plot_micro_evolution_overview(micro_daily, plot_dir)
    plot_micro_genotype_stream(micro_genotype, plot_dir)
    plot_micro_hospital_heatmap(micro_hosp, plot_dir)
    plot_micro_patient_phase_space(micro_patient, plot_dir)

    if quiet:
        print(f"plots_written to={plot_dir}")
    else:
        _print_sanity_checks(daily, hosp, micro_daily)
        if interactive:
            print("Opening interactive micro explorer...")
            launch_micro_patient_explorer(micro_patient)


def main() -> None:
    parser = argparse.ArgumentParser(description="Visualize simulation results.")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--plot-dir", type=Path, default=DEFAULT_PLOT_DIR)
    parser.add_argument(
        "--interactive",
        action="store_true",
        help="Open an interactive day-slider explorer for micro patient phase space.",
    )
    args = parser.parse_args()
    run(data_dir=args.data_dir, plot_dir=args.plot_dir, interactive=args.interactive)


if __name__ == "__main__":
    main()
