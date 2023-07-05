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
        A,B = self.topology.split(self.bond_a,self.indexes) #TODO make sure indexes is being used
        bond_b_in_B = B.get_bond(self.bond_b.atom_a.atom_name, self.bond_b.atom_b.atom_name)
        
        if not bond_b_in_B: # B is on the LHS so swap A and B
            B,A = A,B
            pseudoatom_a = A.pseudoatoms[self.indexes[0]]    #TODO make sure indexes is being used
            pseudoatom_a.index = self.indexes[1]
            pseudoatom_b = B.pseudoatoms[self.indexes[1]]
            pseudoatom_b.index = self.indexes[0]
            bond_b_in_B = B.get_bond(self.bond_b.atom_a.atom_name, self.bond_b.atom_b.atom_name)
        B,C = B.split(self.bond_b,self.indexes)
        if len(C.pseudoatoms) >1: # B is on the RHS so swap B and C
            B,C = C,B
            pseudoatom_b_left = B.pseudoatoms[self.indexes[1]]
            pseudoatom_b_left.index = self.indexes[0]
            pseudoatom_b_right = B.pseudoatoms[self.indexes[0]]
            pseudoatom_b_right.index = self.indexes[1]
            pseudoatom_c = C.pseudoatoms[self.indexes[0]]
            pseudoatom_c.index = self.indexes[0]
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
