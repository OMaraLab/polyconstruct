from __future__ import annotations
import copy
import json
from typing import Tuple, Union, List

from .Junction import Junction, Junctions

from .Bonds import Bond
from .Topology import Topology
    
class Monomer:
    """
    Represents a monomer with its properties in a molecular system.
    """
    def __init__(self, topology: Topology, junctions: Union[Junctions, List[Junction]]):
        """
        Represents a monomer with its properties in a molecular system.

        :param topology: the Topology object used to create this Monomer.
        :type topology: Topology
        :param junctions: a list of Junction objects or Junctions object
                containing all Junctions which will be present in this Monomer.
        :type junctions: Union[Junctions, List[Junction]]
        """
        self.topology = topology
        if isinstance(junctions, Junctions):
            self.junctions = junctions
        else:
            junctions_obj = Junctions()
            for junction in junctions:
                junctions_obj.add(junction)
            self.junctions = junctions_obj

    from .Polymer import Polymer

    @classmethod
    def from_Polymer(cls, polymer: Polymer) -> Monomer:
        """
        Convert a Polymer object to a Monomer object. Useful for creating
        Polymeric branches of a set size and converting them back to a Monomer
        so that they can be polymerised onto another Polymer.

        :param polymer: a Polymer to convert.
        :type polymer: Polymer
        :raises ValueError: if the monomer atom in one of the Polymer's
                Junctions is not present.
        :raises ValueError: if the residue atom in one of the Polymer's
                Junctions is not present.
        :return: a new Monomer which has the same attributes
                (Topology and Junctions) as the provided Polymer.
        :rtype: Monomer
        """
        topology_copy = polymer.topology.copy()
        junctions_copy = []
        for junction in polymer.junctions:
            monomer_atom_id = junction.monomer_atom.atom_id
            residue_atom_id = junction.residue_atom.atom_id
            monomer_atom = topology_copy.get_atom(monomer_atom_id)
            residue_atom = topology_copy.get_atom(residue_atom_id)
            if monomer_atom is None:
                raise ValueError(f"Could not find atom {monomer_atom_id} for junction")
            if residue_atom is None:
                raise ValueError(f"Could not find atom {residue_atom_id} for junction")
            junctions_copy.append(Junction(monomer_atom, residue_atom, junction.name))
        return cls(topology_copy, junctions_copy)

    def copy(self) -> Monomer:
        """
        Method to duplicate a Monomer, which leverages the Topology.copy() function.

        :raises ValueError: if the monomer atom in one of the Polymer's
                Junctions is not present.
        :raises ValueError: if the residue atom in one of the Polymer's
                Junctions is not present.
        :return: a new Monomer which has the same attributes as this Monomer.
        :rtype: Monomer
        """
        new_topology = self.topology.copy()
        new_junctions = Junctions()
        for junction in self.junctions:
            monomer_atom_id = junction.monomer_atom.atom_id
            residue_atom_id = junction.residue_atom.atom_id
            monomer_atom = new_topology.get_atom(monomer_atom_id)
            residue_atom = new_topology.get_atom(residue_atom_id)
            if monomer_atom is None:
                raise ValueError(f"Could not find atom {monomer_atom_id}")
            if residue_atom is None:
                raise ValueError(f"Could not find atom {residue_atom_id}")
            new_junctions.add(Junction(monomer_atom, residue_atom, junction.name))
        new_monomer = Monomer(new_topology, new_junctions)
        return new_monomer

    def renumber_atoms(self, start: int):
        """
        Renumber the ids of Atoms in the Topology starting from a given number.
        Each atom's new id number is equal to its current id number plus the
        value of start. 
        
        E.g. an atom with and id of 1 renumbered with a start value of 10 will
            have a new id of 11.

        :param start: value that will be added to the atom's existing ids to
                renumber them.
        :type start: int
        """
        self.topology.renumber_atoms(start)

    def to_dict(self) -> dict:
        """
        Convert this Monomer to a dictionary representation.

        The structure of the dictionary is as below:
        {"topology": self.topology.to_dict(),
        "junctions": self.junctions.to_dict()}

        :return: a dictionary containing references to the dictionary
                representations of this Monomer's Topology and Junctions
                attributes.
        :rtype: dict
        """
        return{
            "topology": self.topology.to_dict(),
            "junctions": self.junctions.to_dict()
        }

    @classmethod
    def from_dict(cls, data: dict) -> Monomer:
        """
        Create a new Monomer from a dictionary, such as that created with
        Monomer.to_dict().

        The structure of the dictionary is as below:
        {"topology": self.topology.to_dict(),
        "junctions": self.junctions.to_dict()}

        :param data: dictionary containing data to make a Monomer, generate
                with 'to_dict()'.
        :type data: dict
        :return: a new Monomer
        :rtype: Monomer
        """
        topology = Topology.from_dict(data["topology"])
        atoms = topology.atoms
        junctions = Junctions.from_dict(data["junctions"], atoms)
        return cls(topology, junctions)

    def save(self, file_path: str):
        """
        Save and export the Monomer to a JSON text dump.

        :param file_path: path to and desired name of output file.
        :type file_path: str
        """
        with open(file_path, "w") as f:
            json.dump(self.to_dict(), f)

    @classmethod
    def load(cls, file_path: str) -> Monomer:
        """
        Load a JSON text dump of a Monomer, generated with Monomer.save(), to
        a new Monomer.

        :param file_path: path to the Monomer JSON dictionary text dump file.
        :type file_path: str
        :return: a new Monomer
        :rtype: Monomer
        """
        with open(file_path, "r") as f:
            data = json.load(f)

        return cls.from_dict(data)
