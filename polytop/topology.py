import copy
from functools import singledispatchmethod
import re
import warnings
from typing import List, Optional, Tuple

import numpy as np
import numpy.linalg as la

from .angles import Angle
from .atoms import Atom
from .bonds import Bond
from .dihedrals import Dihedral, Dihedral_type
from .exclusions import Exclusion
from .molecule_type import MoleculeType
from .pairs import Pair


class Topology:
    """
    Represents the topology of a molecular system, including atoms, bonds, angles, and dihedrals.

    Parameters
    ----------
    atoms : List[Atom], optional
        A list of atoms in the molecular system. Default is an empty list.
    preamble : str, optional
        A string containing the preamble of the topology file. Default is an empty string.
    molecule_type : str, optional
        The type of molecule represented by the topology. Default is an empty string.

    Attributes
    ----------
    atoms : List[Atom]
        A list of atoms in the molecular system.
    bonds : List[Bond]
        A list of bonds in the molecular system.
    angles : List[Angle]
        A list of angles in the molecular system.
    dihedrals : List[Dihedral]
        A list of dihedral angles in the molecular system.

    Methods
    -------
    from_ITP(cls, itp_file: str) -> 'Topology'
        Class method to create a Topology object from a GROMACS ITP file.

    """
    def __init__(
        self,
        atoms: Optional[List[Atom]] = None,
        preamble: Optional[List[str]] = None,
        molecule_type: Optional[MoleculeType] = None,
    ):
        self.molecule_type = molecule_type
        self.atoms = atoms or []
        self.preamble = preamble or []
        self.reorder_atoms()
        
    def copy(self):
        new_topology = Topology(
            atoms=copy.deepcopy(self.atoms),
            preamble=copy.deepcopy(self.preamble),
            molecule_type=copy.deepcopy(self.molecule_type)
        )
        # add bonds, angles, dihedrals, pairs and exclusions
        for bond in self.bonds:
            new_topology.bonds.append(copy.deepcopy(bond))
        for angle in self.angles:
            new_topology.angles.append(copy.deepcopy(angle))
        for dihedral in self.dihedrals:
            new_topology.dihedrals.append(copy.deepcopy(dihedral))
        for pair in self.pairs:
            new_topology.pairs.append(copy.deepcopy(pair))
        for exclusion in self.exclusions:
            new_topology.exclusions.append(copy.deepcopy(exclusion))
        return new_topology
            
    @property
    def name(self) -> str:
        return "Unknown" if self.molecule_type is None else self.molecule_type.name

    def __repr__(self) -> str:
        return f"Topology: {self.name} ({len(self.atoms)} atoms)"

    @property
    def pseudoatoms(self) -> List[Atom]:
        return {atom.index: atom for atom in self.atoms if atom.is_virtual}

    @property
    def bonds(self) -> List[Bond]:
        bonds = []
        for atom in self.atoms:
            for bond in atom.bonds:
                if bond not in bonds:
                    bonds.append(bond)
        return bonds

    @property
    def angles(self) -> List[Angle]:
        angles = []
        for atom in self.atoms:
            for bond in atom.bonds:
                for angle in bond.angles:
                    if angle not in angles:
                        angles.append(angle)
        return angles

    @property
    def dihedrals(self) -> List[Dihedral]:
        dihedrals = []
        for atom in self.atoms:
            for bond in atom.bonds:
                for angle in bond.angles:
                    for dihedral in angle.dihedrals:
                        if dihedral not in dihedrals:
                            dihedrals.append(dihedral)
        return dihedrals

    @property
    def pairs(self) -> List[Pair]:
        pairs = []
        for atom in self.atoms:
            for pair in atom.pairs:
                if pair not in pairs:
                    pairs.append(pair)
        return pairs

    @property
    def exclusions(self) -> List[Exclusion]:
        exclusions = []
        for atom in self.atoms:
            for exclusion in atom.exclusions:
                if exclusion not in exclusions:
                    exclusions.append(exclusion)
        return exclusions

    @classmethod
    def from_ITP(cls, file_path: str)->'Topology':
        """
        Class method to create a Topology object from a GROMACS ITP file.

        Parameters
        ----------
        itp_file : str
            The path to the GROMACS ITP file.

        Returns
        -------
        topology : Topology
            The created Topology object.

        """        
        with open(file_path, "r") as f:
            lines = f.readlines()

        preamble = []
        molecule_type = None
        atoms = []

        section = None
        for line in lines:
            line = line.strip()
            if not line or line.startswith(";"):
                if section is None:
                    preamble.append(line)
                continue
            if line.startswith("["):
                section = line.strip("[] ").lower()
                continue
            if section == "moleculetype":
                molecule_type = MoleculeType.from_line(line)
                continue
            elif section == "atoms":
                atoms.append(Atom.from_line(line))
            elif section == "bonds":
                Bond.from_line(line, atoms)
            elif section == "angles":
                Angle.from_line(line, atoms)  ##
            elif section == "pairs":
                Pair.from_line(line, atoms)
            elif section == "dihedrals":
                Dihedral.from_line(line, atoms)
            elif section == "exclusions":
                Exclusion.from_line(line, atoms)
            else:
                warnings.warn(f"Unknown section {section} in {file_path}")
        return cls(atoms, preamble, molecule_type)

    def to_ITP(self, file_path: str):
        """
        Write the topology to a GROMACS ITP file.
        Args:
            file_path (str): The path to the GROMACS ITP file.
        """
        with open(file_path, "w") as f:
            for line in self.preamble:
                f.write(line + "\n")

            f.write("\n[ moleculetype ]\n")
            f.write(str(self.molecule_type) + "\n")

            f.write("[ atoms ]\n")
            for atom in self.atoms:
                f.write(str(atom) + "\n")

            f.write("\n[ bonds ]\n")
            for bond in self.bonds:
                f.write(str(bond) + "\n")

            f.write("\n[ pairs ]\n")
            for pair in self.pairs:
                f.write(str(pair) + "\n")

            f.write("\n[ angles ]\n")
            for angle in self.angles:
                f.write(str(angle) + "\n")

            if any(
                dihedral.dihedral_type is Dihedral_type.improper.value
                for dihedral in self.dihedrals
            ):
                f.write("\n[ dihedrals ]\n")
                f.write("; GROMOS improper dihedrals\n")
                for dihedral in self.dihedrals:
                    if dihedral.dihedral_type is Dihedral_type.improper.value:
                        f.write(str(dihedral) + "\n")

            f.write("\n[ dihedrals ]\n")
            for dihedral in self.dihedrals:
                if dihedral.dihedral_type is Dihedral_type.proper.value:
                    f.write(str(dihedral) + "\n")

            f.write("\n[ exclusions ]\n")
            for exclusion in self.exclusions:
                f.write(str(exclusion) + "\n")

    @property
    def netcharge(self):
        # Implementation code here
        return sum(atom.partial_charge for atom in self.atoms)

    @netcharge.setter
    def netcharge(self, new_netcharge):
        """
        NB: The partial charges of all atoms will be adjusted by the same amount.
        To change the partial charge of each atom proportionally (ie: to change the net charge without changing 
        the charge distribution), use the function proportional_charge_change instead.
        """
        old_netcharge = self.netcharge
        charge_difference_per_atom = (new_netcharge - old_netcharge) / len(self.atoms)
        for atom in self.atoms:
            atom.partial_charge += charge_difference_per_atom

    def proportional_charge_change(self, new_netcharge):
        """
        Change the net charge of the molecular system by changing the partial charge of each atom proportionally.
        """
        old_netcharge = self.netcharge
        charge_difference = new_netcharge - old_netcharge
        total_absolute_charge = sum(abs(atom.partial_charge) for atom in self.atoms)
        
        for atom in self.atoms:
            atom.partial_charge += charge_difference * (abs(atom.partial_charge) / total_absolute_charge)
        
    @singledispatchmethod
    def get_atom(self, atom_id: int) -> Atom:
        return next((atom for atom in self.atoms if atom.atom_id == atom_id), None)
    
    @get_atom.register
    def _(self, atom_name: str) -> Atom:
        return next((atom for atom in self.atoms if atom.atom_name == atom_name), None)

    @singledispatchmethod
    def get_bond(self, atom_a_id: int, atom_b_id: int) -> Bond:
        atom_a = self.get_atom(atom_a_id)
        atom_b = self.get_atom(atom_b_id)
        return Bond.from_atoms(atom_a, atom_b)

    @get_bond.register
    def _(self, atom_a_name: str, atom_b_name: str) -> Bond:
        atom_a = self.get_atom(atom_a_name)
        atom_b = self.get_atom(atom_b_name)
        return Bond.from_atoms(atom_a, atom_b)

    @singledispatchmethod
    def get_angle(self, atom_a_id: int, atom_b_id: int, atom_c_id: int) -> Angle:
        atom_a = self.get_atom(atom_a_id)
        atom_b = self.get_atom(atom_b_id)
        atom_c = self.get_atom(atom_c_id)
        return Angle.from_atoms(atom_a, atom_b, atom_c)

    @get_angle.register
    def _(self, atom_a_name: str, atom_b_name: str, atom_c_name: str) -> Angle:
        atom_a = self.get_atom(atom_a_name)
        atom_b = self.get_atom(atom_b_name)
        atom_c = self.get_atom(atom_c_name)
        return Angle.from_atoms(atom_a, atom_b, atom_c)

    @singledispatchmethod
    def get_dihedral(
        self, atom_a_id: int, atom_b_id: int, atom_c_id: int, atom_d_id: int
    ) -> Dihedral:
        atom_a = self.get_atom(atom_a_id)
        atom_b = self.get_atom(atom_b_id)
        atom_c = self.get_atom(atom_c_id)
        atom_d = self.get_atom(atom_d_id)
        return Dihedral.from_atoms(atom_a, atom_b, atom_c, atom_d)

    @get_dihedral.register
    def _(self, atom_a_name: str, atom_b_name: str, atom_c_name: str, atom_d_name: str) -> Dihedral:
        atom_a = self.get_atom(atom_a_name)
        atom_b = self.get_atom(atom_b_name)
        atom_c = self.get_atom(atom_c_name)
        atom_d = self.get_atom(atom_d_name)
        return Dihedral.from_atoms(atom_a, atom_b, atom_c, atom_d)

    def add_atom(self, atom: Atom):
        self.atoms.append(atom)

    def to_dict(self):
        return {
            "atoms": [atom.to_dict() for atom in self.atoms],
            "bonds": [bond.to_dict() for bond in self.bonds],
            "angles": [angle.to_dict() for angle in self.angles],
            "dihedrals": [dihedral.to_dict() for dihedral in self.dihedrals],
            "pairs": [pair.to_dict() for pair in self.pairs],
            "exclusions": [exclusion.to_dict() for exclusion in self.exclusions],
            "preamble": self.preamble,
            "moleculetype": self.molecule_type.to_dict(),
        }

    @classmethod
    def from_dict(cls, data):
        topology = cls()

        for atom_data in data["atoms"]:
            atom = Atom.from_dict(atom_data)
            topology.add_atom(atom)

        for bond_data in data["bonds"]:
            try:
                Bond.from_dict(bond_data,topology.atoms)
            except ValueError:
                pass # if we get an error because a bond has lost an atom, just ignore it

        for angle_data in data["angles"]:
            try:
                Angle.from_dict(angle_data,topology.atoms)
            except ValueError:
                pass # if we get an error because an angle has lost an atom, just ignore it

        for dihedral_data in data["dihedrals"]:
            try:
                Dihedral.from_dict(dihedral_data,topology.atoms)
            except ValueError:
                pass # if we get an error because a dihedral has lost an atom, just ignore it

        for pair_data in data["pairs"]:
            try:
                Pair.from_dict(pair_data,topology.atoms)
            except ValueError:
                pass # if we get an error because a pair has lost an atom, just ignore it

        for exclusion_data in data["exclusions"]:
            try:
                Exclusion.from_dict(exclusion_data,topology.atoms)
            except ValueError:
                pass # if we get an error because an exclusion has lost an atom, just ignore it
        
        for preamble_line in data["preamble"]:
            topology.preamble.append(preamble_line)

        topology.molecule_type = MoleculeType.from_dict(data["moleculetype"])

        return topology

    def remove_atom(self, atom: Atom):
        '''Remove an atom and all bonds, angles and dihedrals associated with it - ignore errors if the atom is not present'''
        try:
            atom.remove()
            self.atoms.remove(atom)        
        except ValueError:
            pass
        

    def reorder_atoms(self):
        for i in range(len(self.atoms)):
            self.atoms[i].atom_id = i + 1  # note atom id's are 1 indexed
            
    def renumber_atoms(self, start):
        for atom in self.atoms:
            atom.atom_id = atom.atom_id + start

    def add(self, topology):
        self.atoms.extend(topology.atoms)
        for bond in topology.bonds:
            new_bond = Bond.from_dict(bond.to_dict(), self.atoms)
            self.bonds.append(new_bond)
        for angle in topology.angles:
            new_angle = Angle.from_dict(angle.to_dict(), self.atoms)
            self.angles.append(new_angle)
        for dihedral in topology.dihedrals:
            try:  # if we get an error because a dihedral has lost an atom, just ignore it
                new_dihedral = Dihedral.from_dict(dihedral.to_dict(), self.atoms)
                self.dihedrals.append(new_dihedral)
            except ValueError:
                pass
        for pair in topology.pairs:
            new_pair = Pair.from_dict(pair.to_dict(), self.atoms)
            self.pairs.append(new_pair)
        for exclusion in topology.exclusions:
            new_exclusion = Exclusion.from_dict(exclusion.to_dict(), self.atoms)
            self.exclusions.append(new_exclusion)

    def reverse(self) -> "Topology":
        copied_atoms = copy.deepcopy(self.atoms)
        reversed_atoms = copied_atoms[::-1]
        return Topology(reversed_atoms, self.preamble, self.molecule_type) 
    
    # def split(self, bond: Bond, indexes: Tuple[int,int]) -> Tuple["Topology", "Topology"]:
    #     lhs_atom, rhs_atom = bond.atom_a, bond.atom_b
    #     if lhs_atom.atom_id > rhs_atom.atom_id:
    #         lhs_atom, rhs_atom = rhs_atom, lhs_atom

    #     lhs_atom_name = lhs_atom.atom_name # store the atom names so that we can find them in the fragments
    #     rhs_atom_name = rhs_atom.atom_name            
            
    #     # Create deep copies of the topology
    #     LHS, RHS = copy.deepcopy(self), copy.deepcopy(self)

    #     LHS_bond = LHS.get_bond(lhs_atom_name, rhs_atom_name)
    #     RHS_bond = RHS.get_bond(lhs_atom_name, rhs_atom_name)

    #     # Remove atoms before lhs_atom_idx and after rhs_atom_idx
    #     LHS_atoms_to_remove = list(LHS_bond.RHS())
    #     for atom in LHS_atoms_to_remove:
    #         if atom != LHS_bond.atom_b:
    #             LHS.remove_atom(atom)
    #     # Replace the respective atoms with virtual atoms
    #     LHS_bond.atom_b.virtualize(indexes[1])
        
    #     RHS_atoms_to_remove = list(RHS_bond.LHS())
    #     for atom in RHS_atoms_to_remove:
    #         if atom != RHS_bond.atom_a:
    #             RHS.remove_atom(atom)
    #     RHS_bond.atom_a.virtualize(indexes[0])                
    #     return LHS, RHS

    # def extend_with_topology(self, extension: "Topology", joints: Tuple[str, str] ):
    #     end_virtual = self.get_atom(joints[1])
    #     if not end_virtual:
    #         raise Exception("No virtual atoms in base topology")
    #     start_virtual = extension.get_atom(joints[0])
    #     if not start_virtual:
    #         raise Exception("No virtual atoms in extension topology")
    #     last_atom = end_virtual.bond_neighbours().pop()
    #     next_atom = start_virtual.bond_neighbours().pop()
    #     last_bond = last_atom.bonds.pop()  # bond backwards into the topology
    #     # last_angle = last_bond.angles.pop()  # angle backwards into the topology
    #     next_bond = next_atom.bonds.pop()  # bond forwards into the extension
    #     # next_angle = next_bond.angles.pop()  # angle forward into the extension
    #     last_bond.atom_b = next_atom
    #     next_bond.atom_a = last_atom
    #     next_atom.bonds.add(last_bond)
    #     last_atom.bonds.add(next_bond)
    #     for atom in self.atoms:
    #         if atom == end_virtual:
    #             self.atoms.remove(atom)
    #     for atom in extension.atoms:
    #         if atom != start_virtual:
    #             self.add_atom(atom)

    def contains_bond(self, candidate: Bond) -> bool:
        return any(bond for bond in self.bonds if bond == candidate)
    
    def __repr__(self) -> str:
        return f"Topology: ({len(self.atoms)} atoms, netcharge={self.netcharge})"
