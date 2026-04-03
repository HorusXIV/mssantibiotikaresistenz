# MSS - Antibiotic Resistance Simulation

MSS is a multi-scale simulation framework for modeling antibiotic resistance dynamics in hospital environments. The codebase is organized around a single `src/mss` package so simulation logic, domain objects, and CLI entry points are separated from configuration, documentation, tests, and generated artifacts.

## Structural Assessment

The previous layout worked functionally but caused avoidable friction:

- Application code lived in several top-level folders (`exchange`, `macro_simulation`, `micro_simulation`) with no single source root.
- Runtime scripts, configuration, generated outputs, and domain modules were mixed at the repository root.
- The filename `simulation.py` meant different things in different places, which made discovery harder.
- Documentation described the simulation model well, but it did not serve as a reliable map of the repository itself.

The repository now follows a clearer rule set:

- All maintained Python code lives under `src/mss/`.
- Root-level folders are reserved for configuration, documentation, tests, containers, and generated artifacts.
- Macro, micro, domain, and CLI responsibilities are separated explicitly.

## Repository Layout

```text
MSS/
├── config/
│   ├── simulation_abx.yml
│   └── simulation_realistic.yml
├── containers/
│   └── mss_image.def
├── docs/
│   ├── organizational/
│   │   └── .gitkeep
│   └── system_overview/
│       ├── Flowchart_v0.mmd
│       ├── Flowchart_v1.mmd
│       ├── MindMap.mmd
│       ├── Pruned_MindMap.mmd
│       ├── amr_system_map.gexf
│       ├── amr_system_map_edges.csv
│       ├── amr_system_map_nodes.csv
│       └── build_gephi_graphs.py
├── logs/
│   ├── .gitkeep
│   └── *.out
├── outputs/
│   ├── csv/
│   │   ├── .gitkeep
│   │   └── *.csv
│   └── plots/
│       ├── .gitkeep
│       └── *.png
├── src/
│   └── mss/
│       ├── __init__.py
│       ├── cli/
│       │   ├── __init__.py
│       │   ├── run_coupled_simulation.py
│       │   └── visualize_results.py
│       ├── domain/
│       │   ├── __init__.py
│       │   └── patient.py
│       └── simulation/
│           ├── __init__.py
│           ├── macro/
│           │   ├── __init__.py
│           │   ├── agents.py
│           │   ├── config.py
│           │   ├── grid.py
│           │   └── simulator.py
│           └── micro/
│               ├── __init__.py
│               ├── engine.py
│               ├── genome.py
│               └── simulator.py
├── tests/
│   ├── __init__.py
│   ├── conftest.py
│   ├── test_run_coupled_simulation.py
│   └── integration_tests/
│       ├── __init__.py
│       ├── test_discharge_turnover.py
│       ├── test_grid.py
│       ├── test_macro_patient_integration.py
│       └── test_micro_patient_integration.py
├── .gitignore
├── .gitlab-ci.yml
├── .pre-commit-config.yaml
├── pyproject.toml
├── README.md
├── report.xml
├── slurm_runner.sh
└── uv.lock
```

## Folder Guide

### `src/mss/`

The only source root for maintained application code.

### `src/mss/domain/`

Shared domain objects used across macro and micro layers.

- `patient.py`: canonical patient model, enums, treatment state, and macro/micro exchange contract.

### `src/mss/simulation/macro/`

Hospital-network simulation logic.

- `config.py`: macro-layer configuration dataclass.
- `agents.py`: Mesa agent wrapper used by the hospital grids.
- `grid.py`: hospital department grid and coarse hospital network grid.
- `simulator.py`: macro simulator orchestration, admissions, transfers, transmission, discharge, and micro coupling.

### `src/mss/simulation/micro/`

Within-host bacterial evolution logic.

- `genome.py`: genome representation, resistance traits, and fitness helpers.
- `engine.py`: micro configuration plus the strain-population evolution engine.
- `simulator.py`: batch-processing interface and episode lifecycle management.

### `src/mss/cli/`

Executable entry points that assemble the application from lower-level modules.

- `run_coupled_simulation.py`: loads YAML configuration, runs the coupled macro/micro simulation, and writes CSV outputs.
- `visualize_results.py`: reads generated CSV outputs and writes plots.

### `config/`

Runtime configuration files. Keep these environment- or scenario-specific, not code-specific.

- `simulation_realistic.yml`: main realistic simulation scenario.
- `simulation_abx.yml`: alternative scenario tuned for antibiotic-focused runs.

### `docs/`

Human-facing project documentation, diagrams, and analytical assets.

- `docs/system_overview/`: system maps, Mermaid diagrams, and graph-building helper script.
- `docs/organizational/`: reserved for planning or process documentation.

### `tests/`

Automated verification for the new `src` layout.

- `conftest.py`: ensures `src/` is on the import path during test runs.
- `test_run_coupled_simulation.py`: configuration loading and initial population setup.
- `integration_tests/`: cross-module behavioral tests for macro, micro, and grid interactions.

### `outputs/`

Generated simulation artifacts. These are not source files.

- `outputs/csv/`: simulation result tables.
- `outputs/plots/`: rendered diagnostic plots.

### `logs/`

Execution logs, including Slurm job output.

### `containers/`

Container definitions and related runtime assets.

## Naming and Placement Conventions

Use these rules for all future additions:

- Put all Python application code under `src/mss/`.
- Put cross-cutting business entities in `src/mss/domain/`.
- Put simulation engines under `src/mss/simulation/<layer>/`.
- Put runnable entry points in `src/mss/cli/`.
- Put scenario YAML in `config/`, never beside code modules.
- Put diagrams, architecture notes, and generated maps in `docs/`.
- Put generated data only in `outputs/` or `logs/`, never under `src/` or `tests/`.

Filename conventions:

- Use `config.py` for configuration-only modules.
- Use `simulator.py` for orchestration objects that coordinate a layer.
- Use `engine.py` for computational kernels or lower-level simulation mechanics.
- Use singular names for domain entities such as `patient.py`.
- Prefer descriptive test filenames that mirror the behavior under test.

Import conventions:

- Import from `mss...`, not from relative top-level folders.
- Keep CLI modules thin; move reusable logic into `domain/` or `simulation/`.
- Avoid circular dependencies by making `domain/` independent of `cli/`.

## Running the Project

Install dependencies:

```bash
uv sync
```

Run the coupled simulation:

```bash
uv run mss-run --config config/simulation_realistic.yml
```

Generate plots from existing CSV output:

```bash
uv run mss-visualize --csv-dir outputs/csv --plot-dir outputs/plots
```

Run tests:

```bash
uv run pytest
```

## Guidelines for Adding New Components

When adding new functionality, place it by responsibility rather than by feature name alone.

### Add a new domain entity

- Put it in `src/mss/domain/`.
- Keep it free of CLI concerns and file-system concerns.
- Export it from `src/mss/domain/__init__.py` if it is part of the package-level API.

### Add a new macro behavior

- Configuration fields go in `src/mss/simulation/macro/config.py`.
- Spatial or topology logic goes in `src/mss/simulation/macro/grid.py`.
- Daily orchestration changes go in `src/mss/simulation/macro/simulator.py`.

### Add a new micro behavior

- Genome- or trait-level calculations go in `src/mss/simulation/micro/genome.py`.
- Population evolution logic goes in `src/mss/simulation/micro/engine.py`.
- Batch execution, persistence, or parallelism changes go in `src/mss/simulation/micro/simulator.py`.

### Add a new CLI command

- Create a new module in `src/mss/cli/`.
- Add an entry under `[project.scripts]` in `pyproject.toml`.
- Keep argument parsing local to the CLI module and import reusable logic from lower layers.

### Add a new scenario configuration

- Add a new YAML file under `config/`.
- Name it after the scenario purpose, not after a temporary experiment.
- Reference it from documentation or job runners if it becomes a supported workflow.

### Add a new test

- Unit-style tests go near the top of `tests/`.
- Cross-layer or behavior-flow tests go in `tests/integration_tests/`.
- Match test names to the module or behavior being verified.

## Scalability Guidance

This structure is intended to scale in a controlled way:

- New simulation layers can be added beside `macro/` and `micro/` under `src/mss/simulation/`.
- Additional domain concepts can grow under `src/mss/domain/` without polluting orchestration code.
- CLI commands can grow independently without forcing domain modules to depend on process-level concerns.
- Scenario growth is isolated to `config/` instead of duplicating logic across scripts.

If a subpackage grows beyond roughly 5-7 files, introduce a focused subfolder only when it creates a clearer boundary, for example `src/mss/simulation/macro/policies/`.

## Maintenance Recommendations

To keep the structure healthy over time:

- Reject new top-level code folders unless they are clearly non-source concerns.
- Avoid reintroducing duplicate “runner” logic in multiple files.
- Rename ambiguous modules early if their responsibility broadens.
- Keep generated artifacts out of `src/`, `tests/`, and `docs/`.
- Update this README whenever a new top-level folder, package, or public entry point is added.
- Review imports during code review; `mss...` should remain the default import root.

## Onboarding Recommendations

For new developers, the fastest path is:

1. Read this README for the repository map.
2. Start with `src/mss/domain/patient.py` to understand the shared contract.
3. Read `src/mss/simulation/macro/simulator.py` and `src/mss/simulation/micro/simulator.py` to understand orchestration boundaries.
4. Run `uv run pytest` to validate the environment.
5. Run `uv run mss-run --config config/simulation_realistic.yml` and inspect `outputs/csv/` and `outputs/plots/`.

During onboarding reviews, emphasize one rule: if a file does not belong under `src/mss`, it should only exist at the repository root if it is configuration, documentation, automation, or generated output.
