from __future__ import annotations
from typing import Dict, List, Tuple, Union

import warnings

from .Bonds import Bond
from .Atoms import Atom

class Angle:
    """
    Represents an angle between three atoms in a molecular system.

    :param atom_a: The first atom involved in the angle.
    :type atom_a: Atom
    :param atom_b: The central atom in the angle.
    :type atom_b: Atom
    :param atom_c: The third atom involved in the angle.
    :type atom_c: Atom
    :param angle_type: The type of the angle.
    :type angle_type: int
    :param angle_value: The value of the angle in degrees.
    :type angle_value: float
    :param force_constant: The force constant associated with the angle.
    :type force_constant: float
    :raises ValueError: If there are not bonds present between atoms A and
            B and atoms B and C
    """
    def __init__(self, atom_a: Atom, atom_b: Atom, atom_c: Atom,
            angle_type: int, angle_value: float, force_constant: float) -> None:
        """
        Represents an angle between three atoms in a molecular system.

        :param atom_a: The first atom involved in the angle.
        :type atom_a: Atom
        :param atom_b: The central atom in the angle.
        :type atom_b: Atom
        :param atom_c: The third atom involved in the angle.
        :type atom_c: Atom
        :param angle_type: The type of the angle.
        :type angle_type: int
        :param angle_value: The value of the angle in degrees.
        :type angle_value: float
        :param force_constant: The force constant associated with the angle.
        :type force_constant: float
        :raises ValueError: If there are not bonds present between atoms A and
                B and atoms B and C
        """
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
        print(f"Bond a = {bond_a}")
        print(f"Bond b = {bond_b}")
        print(f"Atom a = {atom_a}, atom b = {atom_b} and atom c = {atom_c}")
        if bond_a is None or bond_b is None:
            return None
        return next((angle for angle in bond_a.angles if angle in bond_b.angles), None)

    def contains_atom(self, atom: Atom) -> bool:
        """
        Check if this Angle contains a given atom.

        :param atom: the Atom you wish to check if it is in this Angle or not
        :type atom: Atom
        :return: True if the Angle contains the given Atom, or False if not.
        :rtype: bool
        """
        return atom in [self.atom_a, self.atom_b, self.atom_c]

    def clone_angle_changing(self, from_atom: Atom, to_atom: Atom) -> Angle:
        """
        Clone the angle, changing the atom that is being replaced. Used during
        the polymer.extend() algorithm to copy and modify angles where a new
        Monomer is joined to the Polymer. 

        :param from_atom: the outgoing Atom, to be replaced
        :type from_atom: Atom
        :param to_atom: the incoming Atom, will replace the position of the
                outgoing Atom in this Angle
        :type to_atom: Atom
        :raises ValueError: if 'from_atom' is not in the Angle
        :return: the new, modified Angle
        :rtype: Angle
        """
        if self.atom_a == from_atom:  # first atom is being replaced
            new_angle = Angle(to_atom, self.atom_b, self.atom_c, self.angle_type, self.angle_value, self.force_constant)
        elif self.atom_b == from_atom:  # second atom is being replaced
            new_angle = Angle(self.atom_a, to_atom, self.atom_c, self.angle_type, self.angle_value, self.force_constant)
        elif self.atom_c == from_atom:  # third atom is being replaced
            new_angle = Angle(self.atom_a, self.atom_b, to_atom, self.angle_type, self.angle_value, self.force_constant)
        else:
            raise ValueError(f"Atom {from_atom} is not in angle {self}")
        return new_angle

    def other_atom(self, atom: Atom) -> List[Atom]:
        """
        Check if the given Atom is in this Angle and return a list of the other
        atoms present in this Angle (i.e. discluding 'atom').

        :param atom: the Atom you wish to check if it is in this Angle or not
        :type atom: Atom
        :raises ValueError: if 'atom' is not in this Angle
        :return: a list of the Atoms in this Angle, not including the Atom
                provided 'atom'. None if 'atom' is not in this Angle.
        :rtype: List[Atom]
        """
        result = [self.atom_a, self.atom_b, self.atom_c]
        if atom not in result:
            raise ValueError(f"Atom {atom} is not in angle {self}")
        result.remove(atom)
        return result

    def remove(self):
        """
        Delete self from all related Dihedrals and Bonds. Used to cleanup and
        remove attributes during Polymer.extend().
        """
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
    def from_dict(cls, data, atoms) -> Angle:
        """
        Create a new Angle from a dictionary (such as that created with
        Angle.to_dict()) and list of Atoms. Will retrieve an existing Angle if
        it already exists between these Atoms.

        The structure of the dictionary is as below:
        {'atom_a': self.atom_a.atom_id,
        'atom_b': self.atom_b.atom_id,
        'atom_c': self.atom_c.atom_id,
        'angle_type': self.angle_type,
        'angle_value': self.angle_value,
        'force_constant': self.force_constant}

        :param data: dictionary containing data to make an Angle, generate
                with 'to_dict()'.
        :type data: dict
        :param atoms: list of Atoms. The list may contain more than 3 atoms, as
                long as the id's of the three atoms specified in the data dict
                are present.
        :type atoms: List[Atom]
        :return: a new Angle
        :rtype: Angle
        """
        atom_a = next((atom for atom in atoms if atom.atom_id == data['atom_a']), None)
        atom_b = next((atom for atom in atoms if atom.atom_id == data['atom_b']), None)
        atom_c = next((atom for atom in atoms if atom.atom_id == data['atom_c']), None)
        # check for existing angle
        angle = Angle.from_atoms(atom_a, atom_b, atom_c)
        if angle is not None:
            return angle
        else:
            return cls(atom_a, 
                       atom_b, 
                       atom_c, 
                       angle_type=data['angle_type'], 
                       angle_value=data['angle_value'], 
                       force_constant=data['force_constant'])
