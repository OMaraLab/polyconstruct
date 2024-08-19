import copy
import json
from typing import Tuple, Union, List

from polytop.junction import Junction, Junctions

from .bonds import Bond
from .topology import Topology
    
class Monomer:
    def __init__(self, topology: Topology, junctions: Union[Junctions, List[Junction]]):
        self.topology = topology
        if isinstance(junctions, Junctions):
            self.junctions = junctions
        else:
            junctions_obj = Junctions()
            for junction in junctions:
                junctions_obj.add(junction)
            self.junctions = junctions_obj

    from polytop.polymer import Polymer

    @classmethod
    def from_Polymer(cls, polymer: Polymer):
        return cls(polymer.topology, polymer.junctions)

    def copy(self):
        new_topology = copy.deepcopy(self.topology)
        new_junctions = Junctions()
        for junction in self.junctions:
            monomer_atom_id = junction.monomer_atom.atom_id
            residue_id = junction.residue_atom.atom_id
            monomer_atom = new_topology.get_atom(monomer_atom_id)
            residue_atom = new_topology.get_atom(residue_id)
            if residue_atom is None:
                raise ValueError(f"Could not find atom {monomer_atom_id}")
            if residue_atom is None:
                raise ValueError(f"Could not find atom {residue_id}")
            new_junctions.add(Junction(monomer_atom, residue_atom, junction.name))
        new_monomer = Monomer(new_topology, new_junctions)
        return new_monomer

    def renumber_atoms(self, start):
        """
            Renumber the atoms.ids in the topology starting from a given starting number.
        """
        self.topology.renumber_atoms(start)

    def to_dict(self):
        return{
            "topology": self.topology.to_dict(),
            "junctions": self.junctions.to_dict()
        }

    @classmethod
    def from_dict(cls, data):
        topology = Topology.from_dict(data["topology"])
        atoms = topology.atoms
        junctions = Junctions.from_dict(data["junctions"], atoms)
        return cls(topology, junctions)

    def save(self, file_path):

        with open(file_path, "w") as f:
            json.dump(self.to_dict(), f)

    @classmethod
    def load(cls, file_path):
        with open(file_path, "r") as f:
            data = json.load(f)

        return cls.from_dict(data)
