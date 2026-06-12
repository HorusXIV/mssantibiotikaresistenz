from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from mss.simulation.micro import SimulationConfig, build_micro_config

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CANONICAL_CONFIG = PROJECT_ROOT / "config" / "simulation_realistic.yml"


def micro_config(**overrides: Any) -> SimulationConfig:
    """Build a complete micro config for tests with explicit overrides."""
    raw = yaml.safe_load(CANONICAL_CONFIG.read_text(encoding="utf-8")) or {}
    micro_raw = dict(raw["micro"])
    micro_raw.update(overrides)
    return build_micro_config(micro_raw, source=f"{CANONICAL_CONFIG}:micro")
