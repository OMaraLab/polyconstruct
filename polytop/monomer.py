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
    def __init__(self, topology: Topology, bond_a: Bond, bond_b: Bond, indexes:Tuple[int,int]=None):
        self.topology = topology
        self.bond_a = bond_a
        self.bond_b = bond_b
        self.indexes = indexes
        if self.indexes is None:
            self.indexes = (0,1)
        self.LHS, self.link, self.RHS = self.build_monomer()

    def build_monomer(self) -> Tuple[Topology, Topology, Topology]:
        self.bond_a.order = 0
        self.bond_b.order = 0
        A,B = self.topology.split(self.bond_a,self.indexes)
        bond_b_in_B = B.get_bond(self.bond_b.atom_a.atom_name, self.bond_b.atom_b.atom_name)
        
        if not bond_b_in_B:
            B,A = A,B
            A.pseudoatoms[0].atom_name,B.pseudoatoms[0].atom_name = B.pseudoatoms[0].atom_name,A.pseudoatoms[0].atom_name
            bond_b_in_B = B.get_bond(self.bond_b.atom_a.atom_name, self.bond_b.atom_b.atom_name)
        B,C = B.split(self.bond_b,self.indexes)
        if len(C.pseudoatoms) >1:
            B,C = C,B
            B.pseudoatoms[0].atom_name,C.pseudoatoms[0].atom_name = C.pseudoatoms[0].atom_name,B.pseudoatoms[0].atom_name
        
        return A, B, C


    def to_dict(self):
        return{
            "topology": self.topology.to_dict(),
            "bond_a": {'a':self.bond_a.atom_a.atom_name,'b':self.bond_a.atom_b.atom_name},
            "bond_b": {'a':self.bond_b.atom_a.atom_name,'b':self.bond_b.atom_b.atom_name},
            "indexes": {'a':self.indexes[0],'b':self.indexes[1]}
        }

    @classmethod
    def from_dict(cls, data):
        topology = Topology.from_dict(data["topology"])
        bond_a = topology.get_bond(data["bond_a"]["a"],data["bond_a"]["b"])
        bond_b = topology.get_bond(data["bond_b"]["a"],data["bond_b"]["b"])
        indexes = (data['indexes']['a'],data['indexes']['b'])
        return cls(topology, bond_a, bond_b, indexes)

    def save(self, file_path):

        with open(file_path, "w") as f:
            json.dump(self.to_dict(), f)

    @classmethod
    def load(cls, file_path):
        with open(file_path, "r") as f:
            data = json.load(f)

        return cls.from_dict(data)
