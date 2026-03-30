"""Visualize simulation results from outputs/csv/.

Usage:
    python visualize_results.py
    python visualize_results.py --csv-dir outputs/csv --plot-dir outputs/plots

Produces five diagnostic plots that show whether the simulation is behaving
realistically:
  01_population_overview.png   — population size, prevalence, resistance, clinical load
  02_occupancy_stability.png   — per-hospital occupancy (validates Poisson ↔ discharge)
  03_hospital_prevalence.png   — per-hospital carrier prevalence (spatial distribution)
  04_department_mix.png        — Ward / ICU patient counts + isolated patients over time
  05_hospital_grid_heatmap.png — 3×2 grid coloured by final-day carrier prevalence
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DEFAULT_CSV_DIR = Path(__file__).resolve().parent / "outputs" / "csv"
DEFAULT_PLOT_DIR = Path(__file__).resolve().parent / "outputs" / "plots"

# Expected endemic prevalence band (for reference lines)
ENDEMIC_LOW = 0.15
ENDEMIC_HIGH = 0.25

ROLLING_WINDOW = 7  # days for smoothing


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _load(csv_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    daily_path = csv_dir / "macro_daily.csv"
    hosp_path = csv_dir / "macro_daily_by_hospital.csv"
    if not daily_path.exists():
        raise FileNotFoundError(f"CSV not found: {daily_path}\nRun the simulation first.")
    daily = pd.read_csv(daily_path)
    hosp = pd.read_csv(hosp_path) if hosp_path.exists() else pd.DataFrame()
    return daily, hosp


_quiet_mode = False  # set by run() before plotting


def _save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    if not _quiet_mode:
        print(f"  saved → {path}")


def _rolling(series: pd.Series, window: int = ROLLING_WINDOW) -> pd.Series:
    return series.rolling(window, min_periods=1, center=True).mean()


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
    ax.set_xlabel("Day")
    ax.set_ylabel("Patients")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # --- [0,1] Carrier prevalence ---
    ax = axes[0, 1]
    prev = daily["prevalence"] * 100
    ax.plot(days, prev, color="firebrick", lw=1, alpha=0.5, label="Daily")
    ax.plot(days, _rolling(prev), color="firebrick", lw=2, label=f"{ROLLING_WINDOW}d mean")
    ax.set_title("Carrier Prevalence (%)")
    ax.set_xlabel("Day")
    ax.set_ylabel("Prevalence (%)")
    ax.set_ylim(0, max(prev.max() * 1.15, 30))
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # --- [1,0] Resistance evolution ---
    ax = axes[1, 0]
    res = daily["avg_resistant_fraction"]
    ax.plot(days, res, color="darkorange", lw=1, alpha=0.5, label="Daily")
    ax.plot(days, _rolling(res), color="darkorange", lw=2, label=f"{ROLLING_WINDOW}d mean")
    ax.set_title("Avg Resistant Fraction (carriers only)")
    ax.set_xlabel("Day")
    ax.set_ylabel("Resistant fraction")
    ax.set_ylim(0, 1)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # --- [1,1] Clinical load ---
    ax = axes[1, 1]
    ax.plot(days, daily["isolated_count"], color="purple", lw=1.5, label="Isolated (flag)")
    ax.plot(days, daily["abx_on_count"], color="teal", lw=1.5, label="On ABX")
    ax.set_title("Clinical Load (Isolation flag & ABX)")
    ax.set_xlabel("Day")
    ax.set_ylabel("Patients")
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
    fig.suptitle("Hospital Occupancy Over Time", fontsize=14, fontweight="bold")

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

    ax.set_xlabel("Day")
    ax.set_ylabel("Patients in hospital")
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
    fig.suptitle("Carrier Prevalence per Hospital", fontsize=14, fontweight="bold")

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

    ax.set_xlabel("Day")
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
    fig.suptitle("Patient Distribution by Department + Isolation", fontsize=14, fontweight="bold")

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

    ax.set_xlabel("Day")
    ax.set_ylabel("Patients")
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
        f"Hospital Grid — Carrier Prevalence (day {last_day})", fontsize=13, fontweight="bold"
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
    ax.set_title("Green = low prevalence  |  Red = high prevalence", fontsize=9)

    fig.tight_layout()
    _save(fig, plot_dir / "05_hospital_grid_heatmap.png")


# ---------------------------------------------------------------------------
# Plot 6: Department Grid per Hospital (final day)
# ---------------------------------------------------------------------------


def plot_department_grid(
    hosp: pd.DataFrame,
    plot_dir: Path,
    dept_cols: int = 3,
    dept_rows: int = 2,
    icu_rows: int = 1,
) -> None:
    """Show the 3×2 department zone grid for every hospital on the final day.

    Cells are coloured by patient density (patients per zone cell).
    ICU row = top, Ward rows = bottom, bottom-right cell = Isolation zone.
    """
    if hosp.empty:
        print("  skipped 06_department_grid.png (no per-hospital data)")
        return

    last_day = hosp["day"].max()
    final = hosp[hosp["day"] == last_day].copy()
    hospital_ids = sorted(final["hospital_id"].unique())
    n_hospitals = len(hospital_ids)

    n_icu_cells = dept_cols * icu_rows
    n_ward_cells = dept_cols * (dept_rows - icu_rows)

    ncols_fig = min(n_hospitals, 3)
    nrows_fig = math.ceil(n_hospitals / ncols_fig)
    fig, axes = plt.subplots(nrows_fig, ncols_fig, figsize=(4.5 * ncols_fig, 3.5 * nrows_fig))
    fig.suptitle(
        "Department Grid — Patient Density per Zone (final day)", fontsize=13, fontweight="bold"
    )

    # Normalise axes shape for uniform indexing
    axes = np.array(axes).reshape(nrows_fig, ncols_fig)

    max_density = 1.0
    for hid in hospital_ids:
        row_data = final[final["hospital_id"] == hid]
        if row_data.empty:
            continue
        max_density = max(
            max_density,
            row_data["icu_count"].values[0] / max(n_icu_cells, 1),
            row_data["ward_count"].values[0] / max(n_ward_cells, 1),
        )

    for idx, hid in enumerate(hospital_ids):
        r, c = divmod(idx, ncols_fig)
        ax = axes[r, c]
        row_data = final[final["hospital_id"] == hid]
        if row_data.empty:
            ax.axis("off")
            continue

        ward_ct = int(row_data["ward_count"].values[0])
        icu_ct = int(row_data["icu_count"].values[0])
        ward_density = ward_ct / max(n_ward_cells, 1)
        icu_density = icu_ct / max(n_icu_cells, 1)

        # Build density matrix
        grid = np.zeros((dept_rows, dept_cols))
        labels: list[list[str]] = [[""] * dept_cols for _ in range(dept_rows)]
        for gy in range(dept_rows):
            for gx in range(dept_cols):
                if gy < icu_rows:
                    grid[gy, gx] = icu_density
                    labels[gy][gx] = f"ICU\n~{icu_density:.0f}p"
                else:
                    grid[gy, gx] = ward_density
                    labels[gy][gx] = f"Ward\n~{ward_density:.0f}p"

        im = ax.imshow(grid, cmap="YlOrRd", vmin=0, vmax=max_density, aspect="auto")
        for gy in range(dept_rows):
            for gx in range(dept_cols):
                ax.text(gx, gy, labels[gy][gx], ha="center", va="center", fontsize=9, color="black")

        ax.set_title(hid.replace("hospital_", "H"), fontsize=10)
        ax.set_xticks(range(dept_cols))
        ax.set_yticks(range(dept_rows))
        ax.set_xticklabels([f"Col {i}" for i in range(dept_cols)], fontsize=7)
        ax.set_yticklabels(
            ["ICU" if gy < icu_rows else "Ward" for gy in range(dept_rows)],
            fontsize=7,
        )
        fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04, label="p/cell")

    for idx in range(n_hospitals, nrows_fig * ncols_fig):
        r, c = divmod(idx, ncols_fig)
        axes[r, c].axis("off")

    fig.tight_layout()
    _save(fig, plot_dir / "06_department_grid.png")


# ---------------------------------------------------------------------------
# Sanity-check summary printed to console
# ---------------------------------------------------------------------------


def _print_sanity_checks(daily: pd.DataFrame, hosp: pd.DataFrame) -> None:
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
    print("────────────────────────────────────────────────────────────\n")


# ---------------------------------------------------------------------------
# Core entry point (importable)
# ---------------------------------------------------------------------------


def run(
    csv_dir: Path = DEFAULT_CSV_DIR,
    plot_dir: Path = DEFAULT_PLOT_DIR,
    quiet: bool = False,
) -> None:
    """Generate all diagnostic plots from simulation CSVs.

    Parameters
    ----------
    quiet
        When True (used by run_coupled_simulation.py), suppress all output
        except a single summary line printed at the end.
    """
    global _quiet_mode
    _quiet_mode = quiet

    daily, hosp = _load(csv_dir)

    if not quiet:
        print(f"Reading CSVs from: {csv_dir}")
        print(f"  macro_daily:             {len(daily)} rows")
        print(f"  macro_daily_by_hospital: {len(hosp)} rows")
        print(f"\nGenerating plots → {plot_dir}")

    plot_population_overview(daily, plot_dir)
    plot_occupancy_stability(hosp, plot_dir)
    plot_hospital_prevalence(hosp, plot_dir)
    plot_department_mix(daily, plot_dir)
    plot_hospital_grid_heatmap(hosp, plot_dir)
    plot_department_grid(hosp, plot_dir)

    if quiet:
        print(f"plots_written to={plot_dir}")
    else:
        _print_sanity_checks(daily, hosp)


def main() -> None:
    parser = argparse.ArgumentParser(description="Visualize simulation results.")
    parser.add_argument("--csv-dir", type=Path, default=DEFAULT_CSV_DIR)
    parser.add_argument("--plot-dir", type=Path, default=DEFAULT_PLOT_DIR)
    args = parser.parse_args()
    run(csv_dir=args.csv_dir, plot_dir=args.plot_dir)


if __name__ == "__main__":
    main()
