import re
from typing import Any, List, Optional

from polytop.exclusions import Exclusion
from .bonds import Bond
from .pairs import Pair


class Atom:
    """
    Represents an atom with its properties in a molecular system.

    Attributes
    ----------
    atom_id : int
        The unique identifier of the atom.
    atom_type : str
        The type of the atom, usually based on its element.
    residue_id : int
        The unique identifier of the residue containing the atom.
    residue_name : str
        The name of the residue containing the atom.
    atom_name : str
        The name of the atom, often based on its position within the residue.
    charge_group_num : int
        The charge group the atom belongs to.
    partial_charge : float
        The partial charge of the atom.
    mass : float
        The mass of the atom.

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
    ) -> None:
        self.atom_id = atom_id
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

    @property
    def element(self) -> str:
        return re.sub("[^a-zA-Z]", "", self.atom_name)

    @element.setter
    def element(self, value: str):
        self.atom_name = f"{value}{self.index}"

    @property 
    def index(self) -> int:
        return int(re.sub("[^0-9]", "", self.atom_name))
    
    @index.setter
    def index(self, value: int):
        self.atom_name = f"{self.element}{value}"
    
    @property
    def is_virtual(self):
        return self.atom_type == "X"
    
    def virtualize(self, index: int):
        self.atom_type = "X"
        self.atom_name = f"X{index}"
        self.mass = 0.0
        self.partial_charge = 0.0
    
    @classmethod
    def from_line(cls, line: str) -> "Atom":
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
        return f"{self.atom_id:5} {self.atom_type:5} {self.residue_id:5} {self.residue_name:5} {self.atom_name:5} {self.charge_group_num:5} {self.partial_charge:9.3f} {self.mass:9.4f}"

    def __repr__(self) -> str:
        return f"Atom({self.atom_id}{self.element})->[{','.join([str(bond.other_atom(self).atom_id) for bond in self.bonds])}] charge={self.partial_charge}"

    def remove(self):
        while self.bonds:
            self.bonds.pop().remove()
        while self.pairs:
            self.pairs.pop().remove()
        while self.exclusions:
            self.exclusions.pop().remove()
            
    def bond_neighbor(self, bond: Bond = None):
        if not bond:
            if len(self.bonds) != 1:
                raise ValueError(
                    f"Can't infer a neighbour from an atom with {len(self.bonds)} bonds"
                )
            bond = self.bonds[0]
        return bond.other_atom(self)

    def bond_neighbours(self) -> set["Atom"]:
        return {bond.other_atom(self) for bond in self.bonds}

    def angle_neighbours(self) -> set["Atom"]:
        neighbours = set()
        for bond in self.bonds:
            for angle in bond.angles:
                neighbours.add(angle.atom_a)
                neighbours.add(angle.atom_b)
                neighbours.add(angle.atom_c)
        neighbours.remove(self)
        return neighbours

    def dihedral_neighbours(self) -> set["Atom"]:
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

    def to_dict(self):
        return {
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

    @classmethod
    def from_dict(cls, data):
        return cls(
            atom_id=data["atom_id"],
            atom_name=data["atom_name"],
            atom_type=data["atom_type"],
            residue_id=data["residue_id"],
            residue_name=data["residue_name"],
            mass=data["mass"],
            partial_charge=data["partial_charge"],
            charge_group_num=data["charge_group_num"],
            x = data["x"],
            y = data["y"],
            z = data["z"],            
        )
        
    # def __hash__(self) -> int:
    #     return hash((self.residue_name, self.atom_name))