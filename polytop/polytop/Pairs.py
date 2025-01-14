from __future__ import annotations

from typing import Dict, List, Union

class Atom:
    ...

class Pair:
    """
    Represents interactions between a pair of atoms in a molecular system not
    reflected by bonds.

    :param atom_a: The first atom involved in the pair.
    :type atom_a: Atom
    :param atom_b: The second atom involved in the pair.
    :type atom_b: Atom
    :param pair_type: The type of the pair (e.g., 1-4 interactions).
    :type pair_type: int
    """
    def __init__(self, atom_a: Atom, atom_b: Atom, pair_type: int):
        """
        Represents interactions between a pair of atoms in a molecular system not
        reflected by bonds.

        :param atom_a: The first atom involved in the pair.
        :type atom_a: Atom
        :param atom_b: The second atom involved in the pair.
        :type atom_b: Atom
        :param pair_type: The type of the pair (e.g., 1-4 interactions).
        :type pair_type: int
        """
        self.atom_a = atom_a
        self.atom_b = atom_b
        self.pair_type = pair_type
        atom_a.pairs.add(self)
        atom_b.pairs.add(self)

    @classmethod
    def from_line(cls, line: str, atoms: List[Atom]) -> Pair:
        """
        Class method to construct Pair from the line of an ITP file.

        :param line: the ITP file line
        :type line: str
        :param atoms: list of all Atoms in the Topology 
        :type atoms: List[Atom]
        :return: a new Pair object
        :rtype: Pair
        """
        parts = line.split()
        atom_a = atoms[int(parts[0]) - 1]
        atom_b = atoms[int(parts[1]) - 1]
        pair_type = int(parts[2])
        return cls(atom_a, atom_b, pair_type)

    @staticmethod
    def from_atoms(atom_a: Atom, atom_b: Atom) -> Pair:
        """
        Class method to construct Pair from two Atoms.

        :param atom_a: The first atom involved in the pair.
        :type atom_a: Atom
        :param atom_b: The second atom involved in the pair.
        :type atom_b: Atom
        :return: the new Pair.
        :rtype: Pair
        """
        if atom_a is None or atom_b is None:
            return None
        return next(
            (
                pair
                for pair in atom_a.pairs
                if pair.atom_b == atom_b or pair.atom_a == atom_b
            ),
            None,
        )
    
    def remove(self):
        """
        Delete self from all related Atoms. Used to cleanup and
        remove attributes during Polymer.extend().
        """
        if self in self.atom_a.pairs:
            self.atom_a.pairs.remove(self)
        if self in self.atom_b.pairs:
            self.atom_b.pairs.remove(self)

    def __str__(self):
        return f"{self.atom_a.atom_id:>5} {self.atom_b.atom_id:>5} {self.pair_type:>5}"

    def __repr__(self) -> str:
        return f"Pair({self.atom_a.atom_id}, {self.atom_b.atom_id})"

    def to_dict(self) -> dict:
        """
        Convert this Pair to a dictionary representation.

        The structure of the dictionary is as below:
        {'atom_a': self.atom_a.atom_id,
        'atom_b': self.atom_b.atom_id,
        'pair_type': self.pair_type}

        :return: a dictionary containing the id's of its Atoms and the type of
                this Pair.
        :rtype: dict
        """
        return {
            'atom_a': self.atom_a.atom_id,
            'atom_b': self.atom_b.atom_id,
            'pair_type': self.pair_type,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Union[int, float]], atoms: List[Atom]) -> Pair:
        """
        Create a new Pair from a dictionary (such as that created with
        Pair.to_dict()) and list of Atoms. 

        The structure of the dictionary is as below:
        {'atom_a': self.atom_a.atom_id,
        'atom_b': self.atom_b.atom_id,
        'pair_type': self.pair_type}

        :param data: dictionary containing data to make a Pair, generate
                with 'to_dict()'.
        :type data: dict
        :param atoms: list of Atoms. The list may contain more than 2 atoms, as
                long as the id's of the two atoms specified in the data dict
                are present.
        :type atoms: List[Atom]
        :return: a new Pair
        :rtype: Pair
        """
        atom_a = next((atom for atom in atoms if atom.atom_id == data['atom_a']), None)
        atom_b = next((atom for atom in atoms if atom.atom_id == data['atom_b']), None)
        pair_type = data['pair_type']
        return cls(atom_a, atom_b, pair_type)