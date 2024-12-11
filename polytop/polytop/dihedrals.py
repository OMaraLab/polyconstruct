from __future__ import annotations

import warnings
from enum import IntEnum
from typing import Dict, List, Union

from .angles import Angle


class Atom:
    ...

class Dihedral_type(IntEnum):
    proper = 1  
    improper = 2  
    # add additional dihedral types here
    # remember to update the constraint properties

    @property
    def is_rotational_constraint(self) ->bool:
        '''
        Proper dihedral: constrains torsional rotation around the BC bond
        A -◟B
          /
        C◝- D 
        '''
        return self in [Dihedral_type.proper]

    @property
    def is_planar_constraint(self) ->bool: # The two angles are B-A-C and B-A-D
        '''
        Improper dihedral: constrains orientation of D WRT the CAB plane
            B
            |
        C -◜A◝ - D
        '''
        return self in [Dihedral_type.improper]
    
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
        dihedral_type: Dihedral_type,
        phase_angle: float,
        force_constant: float,
        multiplicity: int,
    ) -> None:
        self.atom_a = atom_a
        self.atom_b = atom_b
        self.atom_c = atom_c
        self.atom_d = atom_d
        self.dihedral_type = Dihedral_type(dihedral_type)
        self.phase_angle = phase_angle
        self.force_constant = force_constant
        self.multiplicity = multiplicity
        if self.dihedral_type.is_rotational_constraint:
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
        elif self.dihedral_type.is_planar_constraint:
            if angle_bac := Angle.from_atoms(atom_b, atom_a, atom_c):
                angle_bac.dihedrals.add(self)
                self.angle_a = angle_bac
            else:
                raise ValueError(f"Could not find angle for dihedral: {self}")
            if angle_bad := Angle.from_atoms(atom_b, atom_a, atom_d):
                angle_bad.dihedrals.add(self)
                self.angle_b = angle_bad
            else:
                raise ValueError(f"Could not find angle for dihedral: {self}")
        else:
            raise ValueError(f"Unknown dihedral type: {dihedral_type}")

    @classmethod
    def from_line(cls, line: str, atoms):
        parts = line.split()
        atom_a = atoms[int(parts[0]) - 1]
        atom_b = atoms[int(parts[1]) - 1]
        atom_c = atoms[int(parts[2]) - 1]
        atom_d = atoms[int(parts[3]) - 1]
        dihedral_type = Dihedral_type(int(parts[4]))
        phase_angle = float(parts[5])
        force_constant = float(parts[6])
        if dihedral_type.is_rotational_constraint:
            multiplicity = int(parts[7])
        elif dihedral_type.is_planar_constraint:
            multiplicity = None # multiplicity is not required for improper dihedrals
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

        angle_bac = Angle.from_atoms(atom_b, atom_a, atom_c)
        angle_bad = Angle.from_atoms(atom_b, atom_a, atom_d)

        if (
            angle_abc and angle_bcd and angle_abc.dihedrals & angle_bcd.dihedrals
        ):  # rotational constraint - proper dihedral
            return angle_abc, angle_bcd
        if (
            angle_bac and angle_bad and angle_bac.dihedrals & angle_bad.dihedrals
        ):  # planar constraint - improper dihedral
            return angle_bac, angle_bad
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
            return None
        common_dihedrals = angle_a.dihedrals & angle_b.dihedrals
        return next(iter(common_dihedrals), None)

    def contains_atom(self, atom: Atom) -> bool:
        return atom in [self.atom_a, self.atom_b, self.atom_c, self.atom_d]
    
    def clone_dihedral_changing(self, from_atom: Atom, to_atom: Atom):
        """ Clone the dihedral, changing the atom that is being replaced """
        if self.atom_a == from_atom:
            new_dihedral = Dihedral(to_atom, self.atom_b, self.atom_c, self.atom_d, self.dihedral_type, self.phase_angle, self.force_constant, self.multiplicity)
        elif self.atom_b == from_atom:
            new_dihedral = Dihedral(self.atom_a, to_atom, self.atom_c, self.atom_d, self.dihedral_type, self.phase_angle, self.force_constant, self.multiplicity)
        elif self.atom_c == from_atom:
            new_dihedral = Dihedral(self.atom_a, self.atom_b, to_atom, self.atom_d, self.dihedral_type, self.phase_angle, self.force_constant, self.multiplicity)
        elif self.atom_d == from_atom:
            new_dihedral = Dihedral(self.atom_a, self.atom_b, self.atom_c, to_atom, self.dihedral_type, self.phase_angle, self.force_constant, self.multiplicity)
        else:
            raise ValueError(f"Atom {from_atom} is not in dihedral {self}")
        return new_dihedral

    def __str__(self):
        if self.dihedral_type.is_rotational_constraint:
            return f"{self.atom_a.atom_id:>5} {self.atom_b.atom_id:>5} {self.atom_c.atom_id:>5} {self.atom_d.atom_id:>5} {self.dihedral_type:>5} {self.phase_angle:>10.4f} {self.force_constant:.4e} {self.multiplicity:>5}"
        elif self.dihedral_type.is_planar_constraint:
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
        # check for existing dihedrals
        dihedral = Dihedral.from_atoms(atom_a, atom_b, atom_c, atom_d)
        if dihedral is not None:
            return dihedral
        else:
            return cls(atom_a, atom_b, atom_c, atom_d, 
                    dihedral_type = data['dihedral_type'], 
                    phase_angle = data['phase_angle'], 
                    force_constant = data['force_constant'], 
                    multiplicity = data['multiplicity'])

