from __future__ import annotations

from typing import Dict, List, Union

class Exclusion:
    """
    Represents non-bonded interactions between two atoms in a molecular system.

    :param atom_a: The first atom involved in the exclusion.
    :type atom_a: Atom
    :param atom_b: The second atom involved in the exclusion.
    :type atom_b: Atom
    """
    def __init__(self, atom_a: "Atom", atom_b: "Atom"):
        """
        Represents exclusions (i.e. non-bonded interactions) between two atoms
        in a molecular system.

        :param atom_a: The first atom involved in the exclusion.
        :type atom_a: Atom
        :param atom_b: The second atom involved in the exclusion.
        :type atom_b: Atom
        """
        self.atom_a = atom_a
        self.atom_b = atom_b
        atom_a.exclusions.add(self)
        atom_b.exclusions.add(self)

    @classmethod
    def from_line(cls, line: str, atoms: List["Atom"], indexes=[0, 1]) -> Exclusion:
        """
        Class method to construct Exclusion from the line of an ITP file and a
        list of all Atom's present in the topology.

        :param line: the ITP file line
        :type line: str
        :param atoms: list of all Atoms in the Topology 
        :type atoms: List["Atom"]
        :param indexes: optional tuple of values to index into elements of the
                line, defaults to [0, 1]
        :type indexes: list, optional
        :return: the new Exclusion
        :rtype: Exclusion
        """
        parts = line.split()
        atom_a = atoms[int(parts[indexes[0]]) - 1]
        atom_b = atoms[int(parts[indexes[1]]) - 1]
        return cls(atom_a, atom_b)

    def remove(self):
        """
        Delete self from all related Atoms. Used to cleanup and
        remove attributes during Polymer.extend().
        """
        if self in self.atom_a.exclusions:
            self.atom_a.exclusions.remove(self)
        if self in self.atom_b.exclusions:
            self.atom_b.exclusions.remove(self)

    def __str__(self):
        return f"{self.atom_a.atom_id:>5} {self.atom_b.atom_id:>5}"

    def __repr__(self) -> str:
        return f"Pair({self.atom_a.atom_id}, {self.atom_b.atom_id})"

    def to_dict(self) -> dict:
        """
        Convert this Exclusion to a dictionary representation.

        The structure of the dictionary is as below:
        {'atom_a': self.atom_a.atom_id,
        'atom_b': self.atom_b.atom_id}

        :return: a dictionary containing the id's of the Exclusion's Atoms.
        :rtype: dict
        """
        return {
            'atom_a': self.atom_a.atom_id,
            'atom_b': self.atom_b.atom_id,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Union[int, float]], atoms: List["Atom"]) -> Exclusion:
        """
        Create a new Exclusion from a dictionary (such as that created with
        Exclusion.to_dict()) and list of Atoms.

        The structure of the dictionary is as below:
        {'atom_a': self.atom_a.atom_id,
        'atom_b': self.atom_b.atom_id}

        :param data: dictionary containing data to make an Exclusion, generate
                with 'to_dict()'.
        :type data: dict
        :param atoms: list of Atoms. The list may contain more than 2 atoms, as
                long as the id's of the two atoms specified in the data dict
                are present.
        :type atoms: List["Atom"]
        :return: a new Exclusion
        :rtype: Exclusion
        """
        atom_a = next((atom for atom in atoms if atom.atom_id == data['atom_a']), None)
        atom_b = next((atom for atom in atoms if atom.atom_id == data['atom_b']), None)
        return cls(atom_a, atom_b)