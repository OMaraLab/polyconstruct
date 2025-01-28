from __future__ import annotations

import warnings
from enum import IntEnum
from typing import Dict, List, Union

from .Angles import Angle

class Dihedral_type(IntEnum):
    """
    Enum to track Dihedral types including proper (1) and improper (2).

    Proper dihedrals also include: Ryckaert-Bellemans (3), Fourier (5),
        proper (multiple) (9), tabulated (8) and restricted (10).
    Improper dihedrals also include: periodic improper (4).

    For more information, see
    `GROMACS documentation <https://manual.gromacs.org/nightly/reference-manual/topologies/topology-file-formats.html#tab-topfile2>`
    and view Table 14.
    """
    proper = 1
    improper = 2
    ryckaert_bellemans = 3
    periodic_improper = 4
    fourier = 5
    multiple = 9
    tabulated = 8
    restricted = 10
    # TODO: add additional dihedral types here
    # but remember to update the constraint properties

    @property
    def is_rotational_constraint(self) -> bool:
        """
        Proper dihedral: constrains torsional rotation around the BC bond

        | A -◟B  |
        |   /    |
        | C◝- D  |

        :return: True if this Dihedral is proper, and False if not.
        :rtype: bool
        """
        return self in [Dihedral_type.proper, Dihedral_type.multiple,
                        Dihedral_type.tabulated, Dihedral_type.restricted]

    @property
    def is_rotational_constraint_with_constants(self) -> bool:
        """
        Proper dihedral: constrains torsional rotation around the BC bond, but
        is defined with 6 constants increase of degrees, energy and multiplicity.

        | A -◟B  |
        |   /    |
        | C◝- D  |

        :return: True if this Dihedral is proper with constants, and False if not.
        :rtype: bool
        """
        return self in [Dihedral_type.ryckaert_bellemans, Dihedral_type.fourier]

    @property
    def is_planar_constraint(self) -> bool:
        """
        Improper dihedral: constrains orientation of D WRT the CAB plane. In
                           other words, the two angles are B-A-C and B-A-D

        |     B      |
        |     |      |
        | C -◜A◝ - D |

        :return: True if this Dihedral is improper, and False if not.
        :rtype: bool
        """
        return self in [Dihedral_type.improper]

    @property
    def is_periodic_planar_constraint(self) -> bool:
        """
        Improper dihedral: constrains orientation of D WRT the CAB plane. In
                           other words, the two angles are B-A-C and B-A-D

        |     B      |
        |     |      |
        | C -◜A◝ - D |

        :return: True if this Dihedral is periodic improper, and False if not.
        :rtype: bool
        """
        return self in [Dihedral_type.periodic_improper]
    
class Dihedral:
    """
    Represents a dihedral angle formed by four atoms in a molecular system.

    :param atom_a: The first atom involved in the dihedral angle.
    :type atom_a: Atom
    :param atom_b: The second atom involved in the dihedral angle.
    :type atom_b: Atom
    :param atom_c: The third atom involved in the dihedral angle.
    :type atom_c: Atom
    :param atom_d: The fourth atom involved in the dihedral angle.
    :type atom_d: Atom
    :param dihedral_type: The type of the dihedral angle (e.g., proper, improper).
    :type dihedral_type: Dihedral_type
    :param phase_angle: The phase angle of the dihedral angle in degrees.
    :type phase_angle: float
    :param force_constant: The force constant associated with the dihedral angle.
    :type force_constant: float
    :param multiplicity: The multiplicity of the dihedral angle.
    :type multiplicity: int
    :raises ValueError: If unable to find an Angle for the Dihedral.
    :raises ValueError: If an unknown Dihedral_type is provided.
    """
    def __init__(
        self,
        atom_a: "Atom",
        atom_b: "Atom",
        atom_c: "Atom",
        atom_d: "Atom",
        dihedral_type: Dihedral_type,
        phase_angle: float = None,
        force_constant: float = None,
        multiplicity: int = None,
        constants: list[float] = None,
        format: str = "gromos",
    ) -> None:
        """
        Represents a dihedral angle formed by four atoms in a molecular system.

        :param atom_a: The first atom involved in the dihedral angle.
        :type atom_a: Atom
        :param atom_b: The second atom involved in the dihedral angle.
        :type atom_b: Atom
        :param atom_c: The third atom involved in the dihedral angle.
        :type atom_c: Atom
        :param atom_d: The fourth atom involved in the dihedral angle.
        :type atom_d: Atom
        :param dihedral_type: The type of the dihedral angle (e.g., proper, improper).
        :type dihedral_type: Dihedral_type
        :param phase_angle: The phase angle of the dihedral angle in degrees.
        :type phase_angle: float
        :param force_constant: The force constant associated with the dihedral angle.
        :type force_constant: float
        :param multiplicity: The multiplicity of the dihedral angle.
        :type multiplicity: int
        :param constants: The constants used to describe the dihedra
                Ryckaert-Bellemans or Fourier potential.
        :type constants: list[int]
        :param format: The forcefield the ITP file is formatted as, options are
                "gromos", "amber", "opls" and "charmm"
        :type format: str, defaults to "gromos" for GROMOS forcefields.
        :raises ValueError: If unable to find an Angle for the Dihedral.
        :raises ValueError: If an unknown Dihedral_type is provided.
        """
        self.atom_a = atom_a
        self.atom_b = atom_b
        self.atom_c = atom_c
        self.atom_d = atom_d
        self.dihedral_type = Dihedral_type(dihedral_type)
        self.phase_angle = phase_angle
        self.force_constant = force_constant
        self.multiplicity = multiplicity
        self.constants = constants
        self.format = format
        # angles and bonds for OPLS forcefield atom order
        if self.format == "opls":
            if self.dihedral_type.is_rotational_constraint or self.dihedral_type.is_rotational_constraint_with_constants:
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
            elif self.dihedral_type.is_planar_constraint or self.dihedral_type.is_periodic_planar_constraint:
                if angle_abc := Angle.from_atoms(atom_a, atom_b, atom_c):
                    angle_abc.dihedrals.add(self)
                    self.angle_a = angle_abc
                else:
                    raise ValueError(f"Could not find angle for dihedral: {self}")
                if angle_abd := Angle.from_atoms(atom_a, atom_b, atom_d):
                    angle_abd.dihedrals.add(self)
                    self.angle_b = angle_abd
                else:
                    raise ValueError(f"Could not find angle for dihedral: {self}")
            else:
                raise ValueError(f"Unknown dihedral type: {dihedral_type}")

        # angles and bonds for GROMOS forcefield atom order
        if self.format == "gromos":
            if self.dihedral_type.is_rotational_constraint or self.dihedral_type.is_rotational_constraint_with_constants:
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
            elif self.dihedral_type.is_planar_constraint or self.dihedral_type.is_periodic_planar_constraint:
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
    def from_line(cls, line: str, atoms, format: str = "gromos") -> Dihedral:
        """
        Class method to construct Dihedral from the line of an ITP file and a
        list of all Atom's present in the topology.

        :param line: the ITP file line
        :type line: str
        :param atoms: list of all Atoms in the Topology 
        :type atoms: List[Atom]
        :return: the new Dihedral
        :rtype: Dihedral
        """
        parts = line.split()
        atom_a = atoms[int(parts[0]) - 1]
        atom_b = atoms[int(parts[1]) - 1]
        atom_c = atoms[int(parts[2]) - 1]
        atom_d = atoms[int(parts[3]) - 1]
        dihedral_type = Dihedral_type(int(parts[4]))
        phase_angle = float(parts[5])
        force_constant = float(parts[6])
        constants = []
        if dihedral_type.is_rotational_constraint or dihedral_type.is_periodic_planar_constraint:
            multiplicity = int(parts[7])
        elif dihedral_type.is_planar_constraint or dihedral_type.is_rotational_constraint_with_constants:
            multiplicity = None # multiplicity is not required for improper dihedrals
        else:
            warnings.warn(f"Unknown dihedral type: {dihedral_type}")

        if dihedral_type.is_rotational_constraint_with_constants:
            constants = [float(parts[5]), float(parts[6]), float(parts[7]), float(parts[8]), float(parts[9])]
            if len(parts) > 10: # only append sixth constant if available (in Rychaert-Bellemans, but not in Fourier dihedral type)
                constants.append(float(parts[10]))
            # set angle and force to None
            phase_angle = None
            force_constant = None

        return cls(
            atom_a,
            atom_b,
            atom_c,
            atom_d,
            dihedral_type,
            phase_angle,
            force_constant,
            multiplicity,
            constants,
            format
        )

    @staticmethod
    def find_angles(
        atom_a: "Atom",
        atom_b: "Atom",
        atom_c: "Atom",
        atom_d: "Atom",
    ) -> tuple[Angle|None, Angle|None]:
        """
        Class method to find Angles present in this Dihedral.

        :param atom_a: The first atom involved in the angle.
        :type atom_a: Atom
        :param atom_b: The central atom in the angle.
        :type atom_b: Atom
        :param atom_c: The third atom involved in the angle.
        :type atom_c: Atom
        

        :param atom_a: The first atom involved in the Dihedral.
        :type atom_a: Atom
        :param atom_b: The second atom involved in the Dihedral.
        :type atom_b: Atom
        :param atom_c: The third atom involved in the Dihedral.
        :type atom_c: Atom
        :param atom_d: The fourth atom involved in the Dihedral.
        :type atom_d: Atom
        :return: a tuple containing the two Angles involved in this Dihedral
        :rtype: tuple[Angle, Angle] or tuple[None, None] if dihedral type is
                neither proper or improper
        """
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

    def contains_atom(self, atom: "Atom") -> bool:
        """
        Check if this Dihedral contains a given Atom.

        :param atom: the Atom you wish to check if it is in this Dihedral or not
        :type atom: Atom
        :return: True if the Dihedral contains the given "Atom", or False if not
        :rtype: bool
        """
        return atom in [self.atom_a, self.atom_b, self.atom_c, self.atom_d]
    
    def other_atoms(self, atom: "Atom") -> List["Atom"]:
        """
        Check if the given Atom is in this Dihedral and return a list of the
        other atoms present in this Dihedral (i.e. discluding 'atom').

        :param atom: the Atom you wish to check if it is in this Dihedral or not
        :type atom: Atom
        :raises ValueError: if 'atom' is not in this Dihedral
        :return: a list of the Atoms in this Dihderal, not including the Atom
                provided 'atom'. None if 'atom' is not in this Angle.
        :rtype: List[Atom]
        """
        result = [self.atom_a, self.atom_b, self.atom_c, self.atom_d]
        if atom not in result:
            raise ValueError(f"Atom {atom} is not in dihedral {self}")
        result.remove(atom)
        return result
    
    def remove(self):
        """
        Delete self from all related Angles. Used to cleanup and
        remove attributes during Polymer.extend().
        """
        if self in self.angle_a.dihedrals:
            self.angle_a.dihedrals.remove(self)
        if self in self.angle_b.dihedrals:
            self.angle_b.dihedrals.remove(self)

    @staticmethod
    def from_atoms(
        atom_a: "Atom",
        atom_b: "Atom",
        atom_c: "Atom",
        atom_d: "Atom",
    ) -> Dihedral:
        """
        Class method to construct Dihedral from four Atoms. There must be at
        least two Angles between these atom pairs that correspond to a valid
        Dihedral_type configuration.

        :param atom_a: The first atom involved in the Dihedral.
        :type atom_a: Atom
        :param atom_b: The second atom involved in the Dihedral.
        :type atom_b: Atom
        :param atom_c: The third atom involved in the Dihedral
        :type atom_c: Atom
        :param atom_d: The fourth atom involved in the Dihedral
        :type atom_d: Atom
        :return: the new Dihedral, or None if the Dihedral_type is
                neither proper or improper
        :rtype: Dihedral
        """
        angle_a, angle_b = Dihedral.find_angles(atom_a, atom_b, atom_c, atom_d)
        if angle_a is None or angle_b is None:
            return None
        common_dihedrals = angle_a.dihedrals & angle_b.dihedrals
        return next(iter(common_dihedrals), None)
    
    def clone_dihedral_changing(self, from_atom: "Atom", to_atom: "Atom") -> Dihedral:
        """
        Clone the dihedral, changing the atom that is being replaced. Used
        during the polymer.extend() algorithm to copy and modify angles where a
        new Monomer is joined to the Polymer.

        :param from_atom: the outgoing "Atom", to be replaced
        :type from_atom: Atom
        :param to_atom: the incoming "Atom", will replace the position of the
                outgoing Atom in this Dihedral
        :type to_atom: Atom
        :raises ValueError: if 'from_atom' is not in the Dihedral
        :return: the new, modified Dihedral
        :rtype: Dihedral
        """
        if self.atom_a == from_atom:
            new_dihedral = Dihedral(to_atom, self.atom_b, self.atom_c, self.atom_d, self.dihedral_type, self.phase_angle, self.force_constant, self.multiplicity, self.constants, self.format)
        elif self.atom_b == from_atom:
            new_dihedral = Dihedral(self.atom_a, to_atom, self.atom_c, self.atom_d, self.dihedral_type, self.phase_angle, self.force_constant, self.multiplicity, self.constants, self.format)
        elif self.atom_c == from_atom:
            new_dihedral = Dihedral(self.atom_a, self.atom_b, to_atom, self.atom_d, self.dihedral_type, self.phase_angle, self.force_constant, self.multiplicity, self.constants, self.format)
        elif self.atom_d == from_atom:
            new_dihedral = Dihedral(self.atom_a, self.atom_b, self.atom_c, to_atom, self.dihedral_type, self.phase_angle, self.force_constant, self.multiplicity, self.constants, self.format)
        else:
            raise ValueError(f"Atom {from_atom} is not in dihedral {self}")
        return new_dihedral

    def __str__(self):
        if self.dihedral_type.is_rotational_constraint or self.dihedral_type.is_periodic_planar_constraint:
            return f"{self.atom_a.atom_id:>5} {self.atom_b.atom_id:>5} {self.atom_c.atom_id:>5} {self.atom_d.atom_id:>5} {self.dihedral_type:>5} {self.phase_angle:>10.4f} {self.force_constant:.4e} {self.multiplicity:>5}"
        elif self.dihedral_type.is_planar_constraint:
            return f"{self.atom_a.atom_id:>5} {self.atom_b.atom_id:>5} {self.atom_c.atom_id:>5} {self.atom_d.atom_id:>5} {self.dihedral_type:>5} {self.phase_angle:>10.4f} {self.force_constant:.4e}"
        elif self.dihedral_type.is_rotational_constraint_with_constants:
            return f"{self.atom_a.atom_id:>5} {self.atom_b.atom_id:>5} {self.atom_c.atom_id:>5} {self.atom_d.atom_id:>5} {self.dihedral_type:>5} {self.constants[0]:>10.3f} {self.constants[1]:>10.3f} {self.constants[2]:>10.3f} {self.constants[3]:>10.3f} {self.constants[4]:>10.3f} {self.constants[5]:>10.3f}"

    def __repr__(self) -> str:
        return f"Dihedral({self.atom_a.atom_id}, {self.atom_b.atom_id}, {self.atom_c.atom_id}, {self.atom_d.atom_id})"

    def to_dict(self) -> dict:
        """
        Convert this Dihedral to a dictionary representation.

        The structure of the dictionary is as below:
        {'atom_a': self.atom_a.atom_id,
        'atom_b': self.atom_b.atom_id,
        'atom_c': self.atom_c.atom_id,
        'atom_d': self.atom_d.atom_id,
        'dihedral_type': self.dihedral_type,
        'phase_angle': self.phase_angle,
        'force_constant': self.force_constant,
        'multiplicity': self.multiplicity}

        :return: a dictionary containing the id's of its Atoms and other
                attributes of this Dihedral.
        :rtype: dict
        """
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
    def from_dict(cls, data: Dict[str, Union[int, float]], atoms: List["Atom"]) -> Dihedral:
        """
        Create a new Dihedral from a dictionary (such as that created with
        Dihedral.to_dict()) and list of Atoms. Will retrieve an existing
        Dihedral if it already exists between these Atoms.

        The structure of the dictionary is as below:
        {'atom_a': self.atom_a.atom_id,
        'atom_b': self.atom_b.atom_id,
        'atom_c': self.atom_c.atom_id,
        'atom_d': self.atom_d.atom_id,
        'dihedral_type': self.dihedral_type,
        'phase_angle': self.phase_angle,
        'force_constant': self.force_constant,
        'multiplicity': self.multiplicity}

        :param data: dictionary containing data to make a Dihedral, generate
                with 'to_dict()'.
        :type data: dict
        :param atoms: list of Atoms. The list may contain more than 4 atoms, as
                long as the id's of the four atoms specified in the data dict
                are present.
        :type atoms: List[Atom]
        :return: a new Dihedral
        :rtype: Dihedral
        """
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

