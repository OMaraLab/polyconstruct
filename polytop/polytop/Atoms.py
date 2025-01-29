from __future__ import annotations
import re
from typing import Any, List, Optional
import warnings

from .Exclusions import Exclusion
from .Bonds import Bond
from .Pairs import Pair


class Atom:
    """
    Represents an atom with its properties in a molecular system.

    :param atom_id: The unique identifier of the atom.
    :type atom_id: int
    :param atom_type: The type of the atom, usually based on its element.
    :type atom_type: str
    :param residue_id: The unique identifier of the residue containing the atom.
    :type residue_id: int
    :param residue_name: The name of the residue containing the atom.
    :type residue_name: str
    :param atom_name: The name of the atom, often based on its position within the residue.
    :type atom_name: str
    :param charge_group_num: The charge group the atom belongs to.
    :type charge_group_num: int
    :param partial_charge: The partial charge of the atom.
    :type partial_charge: float
    :param mass: The mass of the atom.
    :type mass: float
    :param x: The Atom's x position in 3D space, defaults to 0.0
    :type x: float, optional
    :param y: The Atom's y position in 3D space, defaults to 0.0
    :type y: float, optional
    :param z: The Atom's z position in 3D space, defaults to 0.0
    :type z: float, optional
    :param formerly: The atom id of the atom this atom was before
            renumbering, defaults to None
    :type formerly: int, optional
    """
    def __init__(
        self,
        atom_id: int,
        atom_type: str,
        residue_id: int,
        residue_name: str,
        atom_name: str,
        charge_group_num: int,
        partial_charge: float,
        mass: float,
        x: float = 0.0,
        y: float = 0.0,
        z: float = 0.0,
        formerly = None
    ) -> None:
        """
        Represents an atom with its properties in a molecular system.

        :param atom_id: The unique identifier of the atom.
        :type atom_id: int
        :param atom_type: The type of the atom, usually based on its element.
        :type atom_type: str
        :param residue_id: The unique identifier of the residue containing the atom.
        :type residue_id: int
        :param residue_name: The name of the residue containing the atom.
        :type residue_name: str
        :param atom_name: The name of the atom, often based on its position within the residue.
        :type atom_name: str
        :param charge_group_num: The charge group the atom belongs to.
        :type charge_group_num: int
        :param partial_charge: The partial charge of the atom.
        :type partial_charge: float
        :param mass: The mass of the atom.
        :type mass: float
        :param x: The Atom's x position in 3D space, defaults to 0.0
        :type x: float, optional
        :param y: The Atom's y position in 3D space, defaults to 0.0
        :type y: float, optional
        :param z: The Atom's z position in 3D space, defaults to 0.0
        :type z: float, optional
        :param formerly: The atom id of the atom this atom was before
                renumbering, defaults to None
        :type formerly: int, optional
        """
        self.atom_id = atom_id
        self.rank = atom_id
        self.atom_type = atom_type
        self.residue_id = residue_id
        self.residue_name = residue_name
        self.atom_name = atom_name
        self.charge_group_num = charge_group_num
        self.partial_charge = partial_charge
        self.mass = mass
        self.bonds = set()
        self.pairs = set()
        self.exclusions = set()
        self.x = x
        self.y = y
        self.z = z
        self.visited = False
        # when renumbering atoms to extend a polymer we need to keep track of
        # where the atom came from
        self.formerly = formerly

    @property
    def element(self) -> str:
        # compatible with GROMOS 54a7 forcefield, ATB and test files
        element_types = {"H": ["HC", "H", "HS14"], 
                         "O": ["O", "OM", "OA", "OE", "OW", "OEOpt", "OAlc", "OA"], 
                         "C": ["C", "CH0", "CH1", "CH2", "CH3", "CH4", "CH2r", "CR1", "CPos", "CAro"],
                         "N": ["N", "NT", "NL", "NR", "NZ", "NE", "NOpt", "NPri"]}
        element_name = [key for key, val in element_types.items() if self.atom_type in val]
        # if len(element_name) == 0:
        #     warnings.warn(f"Atom type '{self.atom_type}' not supported, attempting to derive element from atom name.")
        #     element_name = self.atom_name[0]
        #     if element_name not in list(element_types.keys()):
        #         warnings.warn(f"Unable to derive element from atom name.")
        # else:
        #     element_name = element_name[0]
        # return element_name

        if len(element_name) == 0:
            if self.atom_name[1:].isnumeric():
                element_name = self.atom_name[0] #TO DO: enable 2 letter elements
                warnings.warn(f"Have extracted element from atom name ({self.atom_name}) instead of type. Element set to {element_name}.")
            else:
                warnings.warn(f"Atom type '{self.atom_type}' not supported, attempting to derive element from atom name.")
                element_name = self.atom_name[0]
                if element_name not in list(element_types.keys()):
                    raise(f"Unable to derive element from atom name.")
                else:
                    warnings.warn(f"Atom element set to {element_name}.")
        else:
            element_name = element_name[0]
        return element_name
    
    #TODO: deconvulute following 3 functions, and atom_name, atom_id and index properties...
    @element.setter
    def element(self, value: str):
        """
        Set the Atom's atom_name attribute to {value}{self.index}. For example,
        a value of 'H' and atom index of 1 will make the atom's name 'H1'.

        :param value: The desired element name for this Atom.
        :type value: str
        """
        self.atom_name = f"{value}{self.index}"

    @property 
    def index(self) -> int:
        """
        Retrieve the atom's index number.

        :return: this Atom's index.
        :rtype: int
        """
        index = re.sub("[^0-9]", "", self.atom_name)
        if index != "":
            return int(index)
        else:
            return self.atom_id
    
    @index.setter
    def index(self, value: int):
        """
        Set the Atom's atom_name attribute to {self.element}{value}. For example,
        an element of 'H' and value 1 will make the atom's name 'H1'.

        :param value: numerical value of the Atom's name - such as its index
                within its elemental type (e.g. H1, H2, H3, etc.).
        :type value: int
        """
        self.atom_name = f"{self.element}{value}"
    
    @property
    def is_virtual(self) -> bool:
        """
        Check if atom is 'virtual' or a dummy.

        :return: True if Atom's atom_type is 'X', but False otherwise
        :rtype: bool
        """
        return self.atom_type == "X"
    
    def virtualize(self, index: int):
        """
        Make atom 'virtual' (i.e. a dummy), by setting its type to X, and
        adjusting the name to match, and its mass and partial charges both to 0.0.

        :param index: desired number for the atom's name and index (e.g. an
                index of 1 will set the Atom's name to 'X1')
        :type index: int
        """
        self.atom_type = "X"
        self.atom_name = f"X{index}"
        self.mass = 0.0
        self.partial_charge = 0.0
    
    @classmethod
    def from_line(cls, line: str) -> "Atom":
        """
        Class method to construct Atom from the line of an ITP file

        :param line: the ITP file line
        :type line: str
        :return: the new Atom
        :rtype: Atom
        """
        parts = line.split()
        atom_id = int(parts[0])
        atom_type = parts[1]
        residue_id = int(parts[2])
        residue_name = parts[3]
        atom_name = parts[4]
        charge_group_num = int(parts[5])
        partial_charge = float(parts[6])
        mass = float(parts[7])
        return cls(
            atom_id,
            atom_type,
            residue_id,
            residue_name,
            atom_name,
            charge_group_num,
            partial_charge,
            mass,
        )

    def __str__(self):
        return f"{self.atom_id:5} {self.atom_type:5} {self.residue_id:5} {self.residue_name:5} {self.atom_name:5} {self.charge_group_num:5} {self.partial_charge:9.6f} {self.mass:9.4f}"

    def __repr__(self) -> str:
        return f"{self.residue_name}.{self.atom_name}->[{','.join([f'{bond.other_atom(self).residue_name}.{bond.other_atom(self).atom_name}' for bond in self.bonds])}]"

    def remove(self):
        """
        Delete self from all related Bonds, Pairs and Exclusions. Used to
        cleanup and remove attributes during Polymer.extend().
        """
        while self.bonds:
            self.bonds.pop().remove()
        while self.pairs:
            self.pairs.pop().remove()
        while self.exclusions:
            self.exclusions.pop().remove()

    def bond_neighbours(self) -> set[Atom]:
        """
        List all the atoms that this atom bonds with.

        :return: set of Atoms that this Atom is bonded to.
        :rtype: set[Atom]
        """
        return {bond.other_atom(self) for bond in self.bonds}

    def deduplicate_bonds(self):
        """
        Remove any bonds from this atom that are duplicates. Used by
        Polymer.extend() to remove duplicated bonds after new ones are made
        during extension. The 'extend' function creates two new identical bonds
        between the existing polymer and incoming monomer, one from the polymer
        to the monomer and the other from the monomer to the polymer, but only
        one of them is needed.
        """
        neighbours = []  # list of all atoms this atom bonds to 
        bonds_to_remove = [] 
        for bond in self.bonds:
            if bond.other_atom(self) in neighbours:
                bonds_to_remove.append(bond)
            else:
                neighbours.append(bond.other_atom(self))
        for bond in bonds_to_remove:
            bond.other_atom(self).bonds.remove(bond)  # remove from the other atom first
            self.bonds.remove(bond)

    def angle_neighbours(self) -> set["Atom"]:
        """
        Find neighbouring Atoms that this Atom produces Angles with.

        :return: set of Atoms that this Atom participates in Angles with.
        :rtype: set[Atom]
        """
        neighbours = set()
        for bond in self.bonds:
            for angle in bond.angles:
                neighbours.add(angle.atom_a)
                neighbours.add(angle.atom_b)
                neighbours.add(angle.atom_c)
        neighbours.remove(self)
        return neighbours

    def dihedral_neighbours(self) -> set["Atom"]:
        """
        Find neighbouring Atoms that this Atom produces Dihedrals with.

        :return: set of Atoms that this Atom participates in Dihedrals with.
        :rtype: set[Atom]
        """
        neighbours = set()
        for bond in self.bonds:
            for angle in bond.angles:
                for dihedral in angle.dihedrals:
                    neighbours.add(dihedral.atom_a)
                    neighbours.add(dihedral.atom_b)
                    neighbours.add(dihedral.atom_c)
                    neighbours.add(dihedral.atom_d)
        if self in neighbours:
            neighbours.remove(self)
        return neighbours

    def to_dict(self) -> dict:
        """
        Convert this Atom to a dictionary representation.

        The structure of the dictionary is as below:
        {"atom_id": self.atom_id,
        "atom_type": self.atom_type,
        "residue_id": self.residue_id,
        "residue_name": self.residue_name,
        "atom_name": self.atom_name,
        "charge_group_num": self.charge_group_num,
        "partial_charge": self.partial_charge,
        "mass": self.mass,
        "x": self.x,
        "y": self.y,
        "z": self.z}
        * Note that the self.formerly attribute will only be included with a
        "formerly" key if the attribute is not None.

        :return: a dictionary containing the attributes of this Atom.
        :rtype: dict
        """
        atom_dict = {
            "atom_id": self.atom_id,
            "atom_type": self.atom_type,
            "residue_id": self.residue_id,
            "residue_name": self.residue_name,
            "atom_name": self.atom_name,
            "charge_group_num": self.charge_group_num,
            "partial_charge": self.partial_charge,
            "mass": self.mass,
            "x": self.x,
            "y": self.y,
            "z": self.z,
        }
        
        # if we have decorated this atom with a formerly attribute, save it
        if hasattr(self, 'formerly'):
            atom_dict["formerly"] = self.formerly
        
        return atom_dict

    @classmethod
    def from_dict(cls, data) -> "Atom":
        """
        Create a new Atom from a dictionary, such as that created with
        Atom.to_dict(). 

        The structure of the dictionary is as below:
        {"atom_id": self.atom_id,
        "atom_type": self.atom_type,
        "residue_id": self.residue_id,
        "residue_name": self.residue_name,
        "atom_name": self.atom_name,
        "charge_group_num": self.charge_group_num,
        "partial_charge": self.partial_charge,
        "mass": self.mass,
        "x": self.x,
        "y": self.y,
        "z": self.z}
        * Note that the "formerly" kew will only be present if it's value
        is not None.

        :param data: dictionary containing data to make an Atom, generate
                with 'to_dict()'.
        :type data: dict
        :return: a new Atom
        :rtype: Atom
        """
        kwargs = {
            "atom_id": data["atom_id"],
            "atom_name": data["atom_name"],
            "atom_type": data["atom_type"],
            "residue_id": data["residue_id"],
            "residue_name": data["residue_name"],
            "mass": data["mass"],
            "partial_charge": data["partial_charge"],
            "charge_group_num": data["charge_group_num"],
            "x": data["x"],
            "y": data["y"],
            "z": data["z"],
        }
        
        if "formerly" in data:  # only if this is part of the dictionary
            kwargs["formerly"] = data["formerly"]
        
        return cls(**kwargs)