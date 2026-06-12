"""Persistent state models for micro patient episodes."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict

import numpy as np

from .engine import StrainPopulation
from .genome import NUM_GENES


@dataclass
class EpisodeState:
    """Persistent state for a patient episode carried between days.

    Attributes:
        episode_id: Stable episode identifier used by the macro patient layer.
        patient_id: Patient identifier associated with the episode.
        population: Current within-host strain population.
        day: Last simulated macro day for the episode.
    """

    episode_id: str
    patient_id: str
    population: StrainPopulation
    day: int = 0

    def to_dict(self) -> Dict[str, Any]:
        """Serialize the episode state for storage.

        Returns:
            JSON-compatible dictionary containing population arrays as lists.
        """
        return {
            "episode_id": self.episode_id,
            "patient_id": self.patient_id,
            "day": self.day,
            "genomes": self.population.genomes.tolist(),
            "populations": self.population.populations.tolist(),
            "lineage_ages": self.population.lineage_ages.tolist(),
            "damage_loads": self.population.damage_loads.tolist(),
            "strain_names": self.population.strain_names,
            "strain_ids": self.population.strain_ids,
            "parent_ids": self.population.parent_ids,
            "donor_ids": self.population.donor_ids,
            "founder_ids": self.population.founder_ids,
            "next_strain_serial": self.population.next_strain_serial,
            "strain_namespace": self.population.strain_namespace,
        }

    def to_payload(self) -> Dict[str, Any]:
        """Serialize the episode state for worker dispatch.

        Returns:
            Dictionary containing copied NumPy arrays and metadata for process
            or thread workers.
        """
        return {
            "episode_id": self.episode_id,
            "patient_id": self.patient_id,
            "day": self.day,
            "genomes": self.population.genomes.copy(),
            "populations": self.population.populations.copy(),
            "lineage_ages": self.population.lineage_ages.copy(),
            "damage_loads": self.population.damage_loads.copy(),
            "strain_names": self.population.strain_names.copy(),
            "strain_ids": self.population.strain_ids.copy(),
            "parent_ids": self.population.parent_ids.copy(),
            "donor_ids": self.population.donor_ids.copy(),
            "founder_ids": self.population.founder_ids.copy(),
            "next_strain_serial": self.population.next_strain_serial,
            "strain_namespace": self.population.strain_namespace,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> EpisodeState:
        """Deserialize an episode state from storage format.

        Args:
            data: Dictionary produced by ``to_dict``.

        Returns:
            Reconstructed episode state.
        """
        genomes = np.array(data["genomes"], dtype=np.float32)
        if genomes.size == 0:
            genomes = genomes.reshape(0, NUM_GENES)
        populations = np.array(data["populations"], dtype=np.float64)
        lineage_ages = np.array(data.get("lineage_ages", []), dtype=np.float64)
        damage_loads = np.array(data.get("damage_loads", []), dtype=np.float64)
        strain_names = [str(name) for name in data.get("strain_names", [])]
        strain_ids = [str(value) for value in data.get("strain_ids", [])]
        parent_ids = [str(value) for value in data.get("parent_ids", [])]
        donor_ids = [str(value) for value in data.get("donor_ids", [])]
        founder_ids = [str(value) for value in data.get("founder_ids", [])]
        return cls(
            episode_id=data["episode_id"],
            patient_id=data["patient_id"],
            population=StrainPopulation(
                genomes=genomes,
                populations=populations,
                lineage_ages=lineage_ages if len(lineage_ages) else None,
                damage_loads=damage_loads if len(damage_loads) else None,
                strain_names=strain_names if strain_names else None,
                strain_ids=strain_ids if strain_ids else None,
                parent_ids=parent_ids if parent_ids else None,
                donor_ids=donor_ids if donor_ids else None,
                founder_ids=founder_ids if founder_ids else None,
                next_strain_serial=int(data.get("next_strain_serial", 0)),
                strain_namespace=str(data.get("strain_namespace", data["episode_id"])),
            ),
            day=data.get("day", 0),
        )

    @classmethod
    def from_payload(cls, data: Dict[str, Any]) -> EpisodeState:
        """Deserialize an episode state from worker payload format.

        Args:
            data: Dictionary produced by ``to_payload``.

        Returns:
            Reconstructed episode state with copied NumPy arrays.
        """
        genomes = np.array(data["genomes"], dtype=np.float32, copy=True)
        if genomes.size == 0:
            genomes = genomes.reshape(0, NUM_GENES)

        return cls(
            episode_id=data["episode_id"],
            patient_id=data["patient_id"],
            population=StrainPopulation(
                genomes=genomes,
                populations=np.array(data["populations"], dtype=np.float64, copy=True),
                lineage_ages=np.array(data["lineage_ages"], dtype=np.float64, copy=True),
                damage_loads=np.array(data["damage_loads"], dtype=np.float64, copy=True),
                strain_names=[str(name) for name in data.get("strain_names", [])] or None,
                strain_ids=[str(value) for value in data.get("strain_ids", [])] or None,
                parent_ids=[str(value) for value in data.get("parent_ids", [])] or None,
                donor_ids=[str(value) for value in data.get("donor_ids", [])] or None,
                founder_ids=[str(value) for value in data.get("founder_ids", [])] or None,
                next_strain_serial=int(data.get("next_strain_serial", 0)),
                strain_namespace=str(data.get("strain_namespace", data["episode_id"])),
            ),
            day=data.get("day", 0),
        )
