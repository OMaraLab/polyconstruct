from __future__ import annotations

from typing import Dict, List, Union


class Atom:
    ...


class Exclusion:
    def __init__(self, atom_a: Atom, atom_b: Atom, exclusion_type: int):
        self.atom_a = atom_a
        self.atom_b = atom_b
        self.exclusion_type = exclusion_type
        atom_a.pairs.add(self)
        atom_b.pairs.add(self)

    @classmethod
    def from_line(cls, line: str, atoms: List[Atom]):
        parts = line.split()
        atom_a = atoms[int(parts[0]) - 1]
        atom_b = atoms[int(parts[1]) - 1]
        exclusion_type = int(parts[2])
        return cls(atom_a, atom_b, exclusion_type)

    def remove(self):
        if self in self.atom_a.exclusions:
            self.atom_a.exclusions.remove(self)
        if self in self.atom_b.exclusions:
            self.atom_b.exclusions.remove(self)

    def __str__(self):
        return f"{self.atom_a.atom_id:>5} {self.atom_b.atom_id:>5} {self.exclusion_type:>5}"

    def __repr__(self) -> str:
        return f"Pair({self.atom_a.atom_id}, {self.atom_b.atom_id})"

    def to_dict(self):
        return {
            'atom_a': self.atom_a.atom_id,
            'atom_b': self.atom_b.atom_id,
            'exclusion_type': self.exclusion_type
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Union[int, float]], atoms: List[Atom]):
        atom_a = next((atom for atom in atoms if atom.atom_id == data['atom_a']), None)
        atom_b = next((atom for atom in atoms if atom.atom_id == data['atom_b']), None)
        exclusion_type = data['exclusion_type']
        return cls(atom_a, atom_b, exclusion_type)