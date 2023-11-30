from __future__ import annotations

from typing import Dict, List, Union


class Atom:
    ...


class Pair:
    """ 
    Represents interactions between a pair of atoms in a molecular system not reflected by bonds.
    Attributes
    ----------
    atom_a : Atom
        The first atom involved in the pair.
    atom_b : Atom
        The second atom involved in the pair.
    pair_type : int
        The type of the pair (e.g., 1-4 interactions).
    """
    def __init__(self, atom_a: Atom, atom_b: Atom, pair_type: int):
        self.atom_a = atom_a
        self.atom_b = atom_b
        self.pair_type = pair_type
        atom_a.pairs.add(self)
        atom_b.pairs.add(self)

    @classmethod
    def from_line(cls, line: str, atoms: List[Atom]):
        parts = line.split()
        atom_a = atoms[int(parts[0]) - 1]
        atom_b = atoms[int(parts[1]) - 1]
        pair_type = int(parts[2])
        return cls(atom_a, atom_b, pair_type)

    def remove(self):
        if self in self.atom_a.pairs:
            self.atom_a.pairs.remove(self)
        if self in self.atom_b.pairs:
            self.atom_b.pairs.remove(self)

    def __str__(self):
        return f"{self.atom_a.atom_id:>5} {self.atom_b.atom_id:>5} {self.pair_type:>5}"

    def __repr__(self) -> str:
        return f"Pair({self.atom_a.atom_id}, {self.atom_b.atom_id})"

    def to_dict(self):
        return {
            'atom_a': self.atom_a.atom_id,
            'atom_b': self.atom_b.atom_id,
            'pair_type': self.pair_type,
        }
    @classmethod
    def from_dict(cls, data: Dict[str, Union[int, float]], atoms: List[Atom]):
        atom_a = next((atom for atom in atoms if atom.atom_id == data['atom_a']), None)
        atom_b = next((atom for atom in atoms if atom.atom_id == data['atom_b']), None)
        pair_type = data['pair_type']
        return cls(atom_a, atom_b, pair_type)