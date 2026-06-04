from __future__ import annotations

import math
from typing import Iterable

from mesa.space import MultiGrid

from mss.domain import Department, HealthState, Patient


class HospitalDepartmentGrid:
    """A per-hospital MultiGrid that groups patients by departments."""

    def __init__(
        self,
        cols: int,
        rows: int,
        icu_rows: int,
        mesa_model,
        rng,
        alpha: float = 0.5,
    ) -> None:
        if cols <= 0 or rows <= 0:
            raise ValueError("Grid dimensions must be positive.")
        if icu_rows < 0 or icu_rows > rows:
            raise ValueError("icu_rows must be between 0 and rows inclusive.")

        self._cols = cols
        self._rows = rows
        self._icu_rows = icu_rows
        self._alpha = alpha
        self._rng = rng
        self._grid = MultiGrid(cols, rows, torus=False)

    def dept_type_at(self, pos: tuple[int, int]) -> Department:
        _, y = pos
        if y < self._icu_rows:
            return Department.ICU
        return Department.WARD

    def dept_positions_for(self, dept: Department) -> list[tuple[int, int]]:
        if dept == Department.ICU:
            row_start, row_end = 0, self._icu_rows
        elif dept == Department.WARD:
            row_start = self._icu_rows
            row_end = self._rows
        else:
            return []

        return [(x, y) for y in range(row_start, row_end) for x in range(self._cols)]

    def assign_to_department(self, agent, dept: Department) -> tuple[int, int]:
        positions = self.dept_positions_for(dept)
        if not positions:
            raise ValueError(f"No positions available for department: {dept}")
        pos = self._rng.choice(positions)
        self._grid.place_agent(agent, pos)
        return pos

    def release(self, agent) -> None:
        if agent.pos is None:
            return
        self._grid.remove_agent(agent)

    def get_all_agents(self) -> list:
        agents = []
        for entry in self._grid.coord_iter():
            if len(entry) == 3:
                cell_content, _, _ = entry
            else:
                cell_content, _ = entry
            agents.extend(cell_content)
        return agents

    def get_all_patients(self) -> list[Patient]:
        return [agent.patient for agent in self.get_all_agents()]

    def cell_stats(self) -> list[dict[str, object]]:
        """Per-cell occupancy for every grid cell (including empty ones).

        Returns one record per cell with its coordinates, department, total
        patient count and carrier count. Enables a true per-cell prevalence
        heatmap (carriers / total), exposing the spatial transmission structure
        that proximity_decay_alpha acts on.
        """
        totals: dict[tuple[int, int], int] = {}
        carriers: dict[tuple[int, int], int] = {}
        for agent in self.get_all_agents():
            pos = agent.pos
            totals[pos] = totals.get(pos, 0) + 1
            if agent.patient.state == HealthState.CARRIER:
                carriers[pos] = carriers.get(pos, 0) + 1

        records: list[dict[str, object]] = []
        for y in range(self._rows):
            for x in range(self._cols):
                total = totals.get((x, y), 0)
                carrier = carriers.get((x, y), 0)
                records.append(
                    {
                        "x": x,
                        "y": y,
                        "department": self.dept_type_at((x, y)).value,
                        "total_patients": total,
                        "carriers": carrier,
                        "prevalence": (carrier / total) if total else 0.0,
                    }
                )
        return records

    def proximity_weight(self, pos_a: tuple[int, int], pos_b: tuple[int, int]) -> float:
        distance = self.chebyshev(pos_a[0], pos_a[1], pos_b[0], pos_b[1])
        return math.exp(-self._alpha * distance)

    @staticmethod
    def chebyshev(x1: int, y1: int, x2: int, y2: int) -> int:
        return max(abs(x1 - x2), abs(y1 - y2))


class HospitalNetworkGrid:
    """Spatial layout of hospitals in a coarse network grid."""

    def __init__(self, hospital_ids: Iterable[str], cols: int = 3, mesa_model=None) -> None:
        hospital_ids = list(hospital_ids)
        if cols <= 0:
            raise ValueError("cols must be positive.")

        self._cols = cols
        self._positions: dict[str, tuple[int, int]] = {}

        for idx, hid in enumerate(hospital_ids):
            x = idx % cols
            y = idx // cols
            self._positions[hid] = (x, y)

    def distance(self, h1: str, h2: str) -> float:
        p1 = self._positions[h1]
        p2 = self._positions[h2]
        return math.hypot(p1[0] - p2[0], p1[1] - p2[1])

    def get_position(self, hid: str) -> tuple[int, int]:
        return self._positions[hid]

    def get_neighbors(self, hid: str) -> list[str]:
        if hid not in self._positions:
            return []
        x0, y0 = self._positions[hid]
        neighbors = []
        for other, (x1, y1) in self._positions.items():
            if other == hid:
                continue
            if abs(x1 - x0) + abs(y1 - y0) == 1:
                neighbors.append(other)
        return neighbors
