from __future__ import annotations

from typing import Dict, List, Union


class Atom:
    ...


class Exclusion:
    def __init__(self, atom_a: Atom, atom_b: Atom):
        self.atom_a = atom_a
        self.atom_b = atom_b
        atom_a.pairs.add(self)
        atom_b.pairs.add(self)

    @classmethod
    def from_line(cls, line: str, atoms: List[Atom], indexes=[0, 1]):
        parts = line.split()
        atom_a = atoms[int(parts[indexes[0]]) - 1]
        atom_b = atoms[int(parts[indexes[1]]) - 1]
        return cls(atom_a, atom_b)

    def remove(self):
        if self in self.atom_a.exclusions:
            self.atom_a.exclusions.remove(self)
        if self in self.atom_b.exclusions:
            self.atom_b.exclusions.remove(self)

    def __str__(self):
        return f"{self.atom_a.atom_id:>5} {self.atom_b.atom_id:>5}"

    def __repr__(self) -> str:
        return f"Pair({self.atom_a.atom_id}, {self.atom_b.atom_id})"

    def to_dict(self):
        return {
            'atom_a': self.atom_a.atom_id,
            'atom_b': self.atom_b.atom_id,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Union[int, float]], atoms: List[Atom]):
        atom_a = next((atom for atom in atoms if atom.atom_id == data['atom_a']), None)
        atom_b = next((atom for atom in atoms if atom.atom_id == data['atom_b']), None)
        return cls(atom_a, atom_b)