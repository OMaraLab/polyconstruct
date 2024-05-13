from __future__ import annotations
from typing import Dict, List, Tuple, Union

import warnings

from polytop.bonds import Bond


class Atom:
    ...


class Angle:
    """
    Represents an angle between three atoms in a molecular system.

    Attributes
    ----------
    atom_a : Atom
        The first atom involved in the angle.
    atom_b : Atom
        The central atom in the angle.
    atom_c : Atom
        The third atom involved in the angle.
    angle_type : str
        The type of the angle.
    angle_value : float
        The value of the angle in degrees.
    force_constant : float
        The force constant associated with the angle.

    """
    def __init__(
        self,
        atom_a: Atom,
        atom_b: Atom,
        atom_c: Atom,
        angle_type: int,
        angle_value: float,
        force_constant: float,
    ) -> None:
        self.atom_a = atom_a
        self.atom_b = atom_b
        self.atom_c = atom_c
        if self.atom_a.atom_id > self.atom_c.atom_id:  # keep atom order ascending
            self.atom_a, self.atom_c = self.atom_c, self.atom_a
        self.angle_type = angle_type
        self.angle_value = angle_value
        self.force_constant = force_constant
        self.bond_ab, self.bond_bc = Angle.find_bonds(atom_a, atom_b, atom_c)
        if self.bond_ab is None or self.bond_bc is None:
            raise ValueError(f"Could not find bonds for angle: {self}")
        self.bond_ab.angles.add(self)
        self.bond_bc.angles.add(self)
        self.dihedrals = set()

    @classmethod
    def from_line(cls, line: str, atoms):
        parts = line.split()
        atom_a = atoms[int(parts[0]) - 1]
        atom_b = atoms[int(parts[1]) - 1]
        atom_c = atoms[int(parts[2]) - 1]
        angle_type = int(parts[3])
        angle_value = float(parts[4])
        force_constant = float(parts[5])
        return cls(atom_a, atom_b, atom_c, angle_type, angle_value, force_constant)

    @staticmethod
    def find_bonds(atom_a: Atom, atom_b: Atom, atom_c: Atom):
        bond_ab = Bond.from_atoms(atom_a, atom_b)
        bond_bc = Bond.from_atoms(atom_b, atom_c)
        return bond_ab, bond_bc

    @staticmethod
    def from_atoms(atom_a: Atom, atom_b: Atom, atom_c: Atom):
        bond_a, bond_b = Angle.find_bonds(atom_a, atom_b, atom_c)
        if bond_a is None or bond_b is None:
            warnings.warn(
                f"Could not find bonds for angle: ({atom_a.atom_id} {atom_b.atom_id} {atom_c.atom_id})"
            )
            return None
        return next((angle for angle in bond_a.angles if angle in bond_b.angles), None)

    def references_atom(self, atom: Atom) -> bool:
        return atom in [self.atom_a, self.atom_b, self.atom_c]

    def other_atom(self, atom: Atom) -> List[Atom]:
        result = [self.atom_a, self.atom_b, self.atom_c]
        if atom not in result:
            raise ValueError(f"Atom {atom} is not in angle {self}")
        result.remove(atom)
        return result

    def remove(self):
        while self.dihedrals:
            self.dihedrals.pop().remove()
        if self in self.bond_ab.angles:
            self.bond_ab.angles.remove(self)
        if self in self.bond_bc.angles:
            self.bond_bc.angles.remove(self)

    def __str__(self):
        return f"{self.atom_a.atom_id:>5} {self.atom_b.atom_id:>5} {self.atom_c.atom_id:>5} {self.angle_type:>5} {self.angle_value:>10.4f} {self.force_constant:.4e}"

    def __repr__(self) -> str:
        return f"Angle({self.atom_a.atom_id}, {self.atom_b.atom_id}, {self.atom_c.atom_id})"

    def to_dict(self):
        data = {
            'atom_a': self.atom_a.atom_id,
            'atom_b': self.atom_b.atom_id,
            'atom_c': self.atom_c.atom_id,
            'angle_type': self.angle_type,
            'angle_value': self.angle_value,
            'force_constant': self.force_constant,
        }
        return data

    @classmethod
    def from_dict(cls, data, atoms):
        atom_a = next((atom for atom in atoms if atom.atom_id == data['atom_a']), None)
        atom_b = next((atom for atom in atoms if atom.atom_id == data['atom_b']), None)
        atom_c = next((atom for atom in atoms if atom.atom_id == data['atom_c']), None)
        angle_type = data['angle_type']
        angle_value = data['angle_value']
        force_constant = data['force_constant']
        return cls(atom_a, atom_b, atom_c, angle_type, angle_value, force_constant)
