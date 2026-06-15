"""Run a long-term convergence simulation (10+ years).

This CLI tool uses config/simulation_realistic_micro.yml as a base and overrides
the simulation duration and population size to observe long-term behavior.
"""

from __future__ import annotations


from mss.cli.run_coupled_simulation import (
    DEFAULT_OUTPUT_DIR,
    PROJECT_ROOT,
    load_coupled_settings,
    run,
)

CONFIG_PATH = PROJECT_ROOT / "config" / "simulation_realistic_micro.yml"


def main() -> None:
    # Load baseline settings
    settings = load_coupled_settings(CONFIG_PATH)

    # Apply overrides for long-term convergence
    settings.run.days = 3650  # 10 years
    settings.run.run_id = "long_term_convergence"

    # Increase population for better statistics
    # Original: 980 + 20 = 1000
    # New: 4900 + 100 = 5000
    settings.population.susceptible_count = 4900
    settings.population.carrier_count = 100

    # Proportionally increase admission rate to maintain turnover
    # Original: 126.0 per 1000 people
    # New: 126.0 * 5 = 630.0
    settings.macro.daily_admission_rate = 630.0

    print("Starting long-term convergence simulation (10 years, 5000 people)...")

    # Run the simulation using the standard runner logic
    output_dir = DEFAULT_OUTPUT_DIR / "LongTermConvergence"
    run(settings, output_dir=output_dir)


if __name__ == "__main__":
    main()
