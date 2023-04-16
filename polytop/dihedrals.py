from __future__ import annotations

import warnings
from enum import IntEnum
from typing import Dict, List, Union

from polytop.angles import Angle


class Atom:
    ...


class Dihedral_type(IntEnum):
    proper = 1  # The two angles are A-B-C and B-C-D
    improper = 2  # The two angles are C-A-D and B-A-C


class Dihedral:
    """
    Represents a dihedral angle formed by four atoms in a molecular system.

    Attributes
    ----------
    atom_a : Atom
        The first atom involved in the dihedral angle.
    atom_b : Atom
        The second atom involved in the dihedral angle.
    atom_c : Atom
        The third atom involved in the dihedral angle.
    atom_d : Atom
        The fourth atom involved in the dihedral angle.
    dihedral_type : Dihedral_type
        The type of the dihedral angle (e.g., proper, improper).
    phase_angle : float
        The phase angle of the dihedral angle in degrees.
    force_constant : float
        The force constant associated with the dihedral angle.
    multiplicity : int
        The multiplicity of the dihedral angle.

    """
    def __init__(
        self,
        atom_a: Atom,
        atom_b: Atom,
        atom_c: Atom,
        atom_d: Atom,
        dihedral_type: int,
        phase_angle: float,
        force_constant: float,
        multiplicity: int,
    ) -> None:
        self.atom_a = atom_a
        self.atom_b = atom_b
        self.atom_c = atom_c
        self.atom_d = atom_d
        self.dihedral_type = dihedral_type
        self.phase_angle = phase_angle
        self.force_constant = force_constant
        self.multiplicity = multiplicity
        if self.dihedral_type == Dihedral_type.proper.value:
            if angle_abc := Angle.from_atoms(atom_a, atom_b, atom_c):
                angle_abc.dihedrals.add(self)
                self.angle_a = angle_abc
            else:
                raise ValueError(f"Could not find angle for dihedral: {self}")
            if angle_bcd := Angle.from_atoms(atom_b, atom_c, atom_d):
                angle_bcd.dihedrals.add(self)
                self.angle_b = angle_bcd
            else:
                raise ValueError(f"Could not find angle for dihedral: {self}")
        elif self.dihedral_type == Dihedral_type.improper.value:
            if angle_cad := Angle.from_atoms(atom_c, atom_a, atom_d):
                angle_cad.dihedrals.add(self)
                self.angle_a = angle_cad
            else:
                raise ValueError(f"Could not find angle for dihedral: {self}")
            if angle_bac := Angle.from_atoms(atom_b, atom_a, atom_c):
                angle_bac.dihedrals.add(self)
                self.angle_b = angle_bac
            else:
                raise ValueError(f"Could not find angle for dihedral: {self}")

    @classmethod
    def from_line(cls, line: str, atoms):
        parts = line.split()
        atom_a = atoms[int(parts[0]) - 1]
        atom_b = atoms[int(parts[1]) - 1]
        atom_c = atoms[int(parts[2]) - 1]
        atom_d = atoms[int(parts[3]) - 1]
        dihedral_type = int(parts[4])
        phase_angle = float(parts[5])
        force_constant = float(parts[6])
        if dihedral_type == Dihedral_type.proper.value:
            multiplicity = int(parts[7])
        elif dihedral_type == Dihedral_type.improper.value:
            multiplicity = None
        else:
            warnings.warn(f"Unknown dihedral type: {dihedral_type}")

        return cls(
            atom_a,
            atom_b,
            atom_c,
            atom_d,
            dihedral_type,
            phase_angle,
            force_constant,
            multiplicity,
        )

    @staticmethod
    def find_angles(
        atom_a: Atom,
        atom_b: Atom,
        atom_c: Atom,
        atom_d: Atom,
    ):
        angle_abc = Angle.from_atoms(atom_a, atom_b, atom_c)
        angle_bcd = Angle.from_atoms(atom_b, atom_c, atom_d)

        angle_cad = Angle.from_atoms(atom_c, atom_a, atom_d)
        angle_bac = Angle.from_atoms(atom_b, atom_a, atom_c)

        if (
            angle_abc and angle_bcd and angle_abc.dihedrals & angle_bcd.dihedrals
        ):  # proper dihedral
            return angle_abc, angle_bcd
        if (
            angle_cad and angle_bac and angle_cad.dihedrals & angle_bac.dihedrals
        ):  # improper dihedral
            return angle_cad, angle_bac
        return None, None

    def references_atom(self, atom: Atom) -> bool:
        return atom in [self.atom_a, self.atom_b, self.atom_c, self.atom_d]
    
    def other_atoms(self, atom: Atom) -> List[Atom]:
        result = [self.atom_a, self.atom_b, self.atom_c, self.atom_d]
        if atom not in result:
            raise ValueError(f"Atom {atom} is not in dihedral {self}")
        result.remove(atom)
        return result
    
    def remove(self):
        if self in self.angle_a.dihedrals:
            self.angle_a.dihedrals.remove(self)
        if self in self.angle_b.dihedrals:
            self.angle_b.dihedrals.remove(self)

    @staticmethod
    def from_atoms(
        atom_a: Atom,
        atom_b: Atom,
        atom_c: Atom,
        atom_d: Atom,
    ):
        angle_a, angle_b = Dihedral.find_angles(atom_a, atom_b, atom_c, atom_d)
        if angle_a is None or angle_b is None:
            warnings.warn(
                f"Could not find angles for dihedral: ({atom_a.atom_id} {atom_b.atom_id} {atom_c.atom_id} {atom_d.atom_id})"
            )
            return None
        common_dihedrals = angle_a.dihedrals & angle_b.dihedrals
        return next(iter(common_dihedrals), None)

    def __str__(self):
        if self.dihedral_type == Dihedral_type.proper:
            return f"{self.atom_a.atom_id:>5} {self.atom_b.atom_id:>5} {self.atom_c.atom_id:>5} {self.atom_d.atom_id:>5} {self.dihedral_type:>5} {self.phase_angle:>10.4f} {self.force_constant:.4e} {self.multiplicity:>5}"
        elif self.dihedral_type == Dihedral_type.improper:
            return f"{self.atom_a.atom_id:>5} {self.atom_b.atom_id:>5} {self.atom_c.atom_id:>5} {self.atom_d.atom_id:>5} {self.dihedral_type:>5} {self.phase_angle:>10.4f} {self.force_constant:.4e}"

    def __repr__(self) -> str:
        return f"Dihedral({self.atom_a.atom_id}, {self.atom_b.atom_id}, {self.atom_c.atom_id}, {self.atom_d.atom_id})"

    def to_dict(self):
        return {
            'atom_a': self.atom_a.atom_id,
            'atom_b': self.atom_b.atom_id,
            'atom_c': self.atom_c.atom_id,
            'atom_d': self.atom_d.atom_id,
            'dihedral_type': self.dihedral_type,
            'phase_angle': self.phase_angle,
            'force_constant': self.force_constant,
            'multiplicity': self.multiplicity,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Union[int, float]], atoms: List[Atom]):
        atom_a = next((atom for atom in atoms if atom.atom_id == data['atom_a']), None)
        atom_b = next((atom for atom in atoms if atom.atom_id == data['atom_b']), None)
        atom_c = next((atom for atom in atoms if atom.atom_id == data['atom_c']), None)
        atom_d = next((atom for atom in atoms if atom.atom_id == data['atom_d']), None)
        dihedral_type = data['dihedral_type']
        phase_angle = data['phase_angle']
        force_constant = data['force_constant']
        multiplicity = data['multiplicity']

        return cls(atom_a, atom_b, atom_c, atom_d, dihedral_type, phase_angle, force_constant, multiplicity)

