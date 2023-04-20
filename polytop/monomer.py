import json
from typing import Tuple

from .bonds import Bond
from .topology import Topology

from enum import Enum
class SubTopologyType(Enum):
    LHS = "LHS"
    LINK = "LINK"
    RHS = "RHS"
    
class Monomer:
    def __init__(self, topology: Topology, bond_a: Bond, bond_b: Bond):
        self.topology = topology
        self.bond_a = bond_a
        self.bond_b = bond_b
        self.LHS, self.link, self.RHS = self.build_monomer()

    def build_monomer(self) -> Tuple[Topology, Topology, Topology]:
        A,B = self.topology.split(self.bond_a,0)
        bond_b_in_B = B.get_bond(self.bond_b.atom_a.atom_name, self.bond_b.atom_b.atom_name)
        
        if not bond_b_in_B:
            B,A = A,B
            bond_b_in_B = B.get_bond(self.bond_b.atom_a.atom_name, self.bond_b.atom_b.atom_name)
        B,C = B.split(self.bond_b,1)
        if len(C.pseudoatoms) >1:
            B,C = C,B
        
        return A, B, C


    def to_dict(self):
        return{
            "topology": self.topology.to_dict(),
            "bond_a": self.bond_a.to_dict(),
            "bond_b": self.bond_b.to_dict(),
        }

    @classmethod
    def from_dict(cls, data):
        topology = Topology.from_dict(data["topology"])
        bond_a = Bond.from_dict(data["bond_a"],topology.atoms)
        bond_b = Bond.from_dict(data["bond_b"],topology.atoms)
        return cls(topology, bond_a, bond_b)

    def save(self, file_path):

        with open(file_path, "w") as f:
            json.dump(self.to_dict(), f)

    @classmethod
    def load(cls, file_path):
        with open(file_path, "r") as f:
            data = json.load(f)

        return cls.from_dict(data)
