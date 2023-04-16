from __future__ import annotations

import json
import random
import re
from typing import Dict, List, Optional, Tuple, Union

from .monomer import Monomer
from .topology import Topology


class Polymer:
    def __init__(
        self,
        monomers: List[Monomer],
        distribution: List[float],
        num_monomers: int,
        seed: Optional[int] = None,
        start_monomer: Optional[Monomer] = None,
        end_monomer: Optional[Monomer] = None,
    ):
        self.monomers = monomers
        self.distribution = distribution
        self.seed = seed
        self.num_monomers = num_monomers
        self.start_monomer = start_monomer
        self.end_monomer = end_monomer

    def get_topology(self) -> Topology:
        polymer_topology = Topology()
        if self.seed:
            random.seed(self.seed)
        else:
            random.seed = random.randint  # unseeded randomization
        # TODO implement build polymer
        monomers_remaining = self.num_monomers

        if self.start_monomer:
            polymer_topology.extend_with_topology(self.start_monomer.LHS)
            monomers_remaining -= 1

        monomers_remaining -= 1

        for _ in range(monomers_remaining):
            chosen_monomer = random.choices(self.monomers, weights=self.distribution)[0]
            # polymer_topology.extend_with_topology(chosen_monomer.link)

        if self.end_monomer:
            polymer_topology.extend_with_topology(self.end_monomer.RHS)
        else:
            chosen_monomer = random.choices(self.monomers, weights=self.distribution)[0]
            polymer_topology.extend_with_topology(chosen_monomer.RHS)

        return polymer_topology

    def save_to_file(self, filename: str) -> None:
        with open(filename, "w") as f:
            json.dump(self.to_dict(), f)

    def load_from_file(self, filename: str) -> None:
        with open(filename, "r") as f:
            data = json.load(f)
        self.from_dict(data)

    def to_dict(
        self,
    ) -> Dict[
        str,
        Union[
            List[Dict[str, Union[float, int]]],
            int,
            Optional[Dict[str, Union[float, int]]],
        ],
    ]:
        return {
            "monomers": [monomer.to_dict() for monomer in self.monomers],
            "distribution": self.distribution,
            "num_monomers": self.num_monomers,
            "seed": self.seed,
            "start_monomer": self.start_monomer.to_dict()
            if self.start_monomer
            else None,
            "end_monomer": self.end_monomer.to_dict() if self.end_monomer else None,
        }

    @classmethod
    def from_dict(
        cls,
        data: Dict[
            str,
            Union[
                List[Dict[str, Union[float, int]]],
                int,
                Optional[Dict[str, Union[float, int]]],
            ],
        ],
    ) -> Polymer:
        monomers = [
            Monomer.from_dict(monomer_data) for monomer_data in data["monomers"]
        ]
        distribution = data["distribution"]
        num_monomers = data["num_monomers"]
        seed = data["seed"]
        start_monomer_data = data.get("start_monomer")
        start_monomer = (
            Monomer.from_dict(start_monomer_data) if start_monomer_data else None
        )
        end_monomer_data = data.get("end_monomer")
        end_monomer = Monomer.from_dict(end_monomer_data) if end_monomer_data else None
        return cls(
            monomers, distribution, seed, num_monomers, start_monomer, end_monomer
        )
