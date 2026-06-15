"""Run an ensemble of simulations to test scenario similarity.

This CLI tool runs the simulation multiple times with different random seeds
to evaluate the variance and stability of the results.
"""

from __future__ import annotations

import copy

from mss.cli.run_coupled_simulation import (
    DEFAULT_OUTPUT_DIR,
    PROJECT_ROOT,
    load_coupled_settings,
    run,
)

CONFIG_PATH = PROJECT_ROOT / "config" / "simulation_realistic_micro.yml"
N_SEEDS = 10


def main() -> None:
    # Load baseline settings
    base_settings = load_coupled_settings(CONFIG_PATH)

    ensemble_dir = DEFAULT_OUTPUT_DIR / "EnsembleConvergence"
    ensemble_dir.mkdir(parents=True, exist_ok=True)

    print(f"Starting ensemble convergence test with {N_SEEDS} seeds...")

    for i in range(1, N_SEEDS + 1):
        print(f"\n--- Run {i}/{N_SEEDS} (seed={i}) ---")

        # Create a deep copy of settings to modify the seed and run_id
        settings = copy.deepcopy(base_settings)
        settings.run.seed = i
        settings.run.run_id = f"ensemble_seed_{i}"
        settings.run.quiet = True  # Keep output clean during ensemble

        # Use a sub-directory for each run
        run_output_dir = ensemble_dir / f"run_seed_{i}"

        run(settings, output_dir=run_output_dir)

    print(f"\nEnsemble test complete. Results saved to: {ensemble_dir}")


if __name__ == "__main__":
    main()
