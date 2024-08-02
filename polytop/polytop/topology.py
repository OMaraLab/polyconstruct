import copy
from functools import singledispatchmethod
import re
import warnings
from typing import Dict, List, Optional, Tuple
# import numpy as np
# import numpy.linalg as la
from .angles import Angle
from .atoms import Atom
from .bonds import Bond
from .junction import Junction, Junctions
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
    from_ITP(cls, itp_file: str, preprocess) -> 'Topology'
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
        self.title = "Unknown molecule"
        if self.preamble and self.preamble[1].startswith(';'):
            self.title = self.preamble[1].lstrip('; ')
        self.reorder_atoms()
        
    def copy(self):
        new_topology = Topology(
            atoms=copy.deepcopy(self.atoms),
            preamble=copy.deepcopy(self.preamble),
            molecule_type=copy.deepcopy(self.molecule_type)
        )
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
    def numerically_order_oxygens(instance: 'Topology', section: str, line: str)->str:
        """
        Preprocesses oxygen atoms with alphanumeric ordering to numeric 
        ordering. For example: OD becomes O4.

        Parameters
        ----------
        instance : Topology
            The Topology class being called.
        section : str
            The itp file section, for example "molecule".
        line : str
            A line from the itp file.

        Returns
        -------
        new_line : str
            The processed line with oxygen atoms renamed numerically.
        """
        if section == 'atoms':
            new_line = re.sub(r'(\bO[A-Z]\b)', lambda match: 'O' + str(ord(match.group(1)[-1]) - ord('A') + 1), line) 
        else:
            new_line = line
        return new_line
    
    @classmethod
    def from_ITP(cls, file_path: str, preprocess=None)->'Topology':
        """
        Class method to create a Topology object from a GROMACS ITP file.

        Parameters
        ----------
        itp_file : str
            The path to the GROMACS ITP file.
        preprocess : lambda function
            Function to preprocess the topology.

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
            if preprocess is not None:
                line = preprocess(section, line)
            if section == "moleculetype":
                molecule_type = MoleculeType.from_line(line)
                continue
            elif section == "atoms":
                atom = Atom.from_line(line)
                # check for issues with atom type or name
                atom.element
                atoms.append(atom)
            elif section == "bonds":
                Bond.from_line(line, atoms)
            elif section == "angles":
                Angle.from_line(line, atoms)  ##
            elif section == "pairs":
                Pair.from_line(line, atoms)
            elif section == "dihedrals":
                Dihedral.from_line(line, atoms)
            elif section == "exclusions":
                if len(line.split()) > 2:
                    for second_atom in range(1, len(line.split())):
                        indexes = [0, second_atom]
                        Exclusion.from_line(line, atoms, indexes=indexes)
                else:
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
    def residue_name(self):
        return self.atoms[0].residue_name
    
    @residue_name.setter
    def residue_name(self, new_name):
        if len(new_name) > 5:
            raise ValueError("Residue name must be 5 characters or less")
        for atom in self.atoms:
            atom.residue_name = new_name
        
    def atom_counts(self) -> Dict[str, int]:
        """ 
        Get the number of atoms of each element in the topology. 
        """
        atom_counts = {}
        for atom in self.atoms:
            if atom.element not in atom_counts:
                atom_counts[atom.element] = 1
            else:
                atom_counts[atom.element] += 1
        return atom_counts
    
    def atom_elements (self) -> set[str]:
        """set of all atom elements in the topology"""
        return set(atom.element for atom in self.atoms)
    
    def auto_rename_atoms(self):
        """
        rename all atoms in the topology to their element symbol plus their index among atoms with the same element (ordered by atom.id)
        
        This compresses the namespace of atom names to the minimum possible while preserving the element and order of atoms in the topology.
        """
        for atom_symbol in self.atom_elements():
            self.reorder_atom_indexes(atom_symbol,1)
        
    def reorder_atom_indexes(self,atom_symbol: str,new_first_index: int):
        """
        Reorder the atom indexes of all atoms with a given symbol to start from a given index.
        """
        atoms_with_this_symbol = [atom for atom in self.atoms if atom.element == atom_symbol]
        atoms_with_this_symbol.sort(key=lambda atom: atom.index)
        for atom in atoms_with_this_symbol:
            atom.index = new_first_index
            atom.element = atom_symbol
            new_first_index += 1

    def max_atom_index(self) -> dict[str, int]:
        """ 
        Get the maximum atom index for each atom symbol in the topology. 
        """
        max_atom_index = {}
        for atom in self.atoms:
            if atom.element not in max_atom_index:
                max_atom_index[atom.element] = atom.index
            else:
                max_atom_index[atom.element] = max(max_atom_index[atom.element], atom.index)
        return max_atom_index
            
    def junction(self, monomer_atom_name:str, residue_atom_name:str, name:str = None) -> Junction:
        monomer_atom = self.get_atom(monomer_atom_name)
        residue_atom = self.get_atom(residue_atom_name)
        return Junction(monomer_atom, residue_atom, name)
    
    @property
    def residue_id(self):
        return self.atoms[0].residue_id
    @residue_id.setter
    def residue_id(self, new_id):
        if new_id < 1:
            raise ValueError("Residue id must be greater than zero")
        if new_id > 99999:
            raise ValueError("Residue id must be less than 100000")
        for atom in self.atoms:
            atom.residue_id = new_id

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
    def get_bond(self, atom_a: Atom, atom_b: Atom) -> Bond:
        return Bond.from_atoms(atom_a, atom_b)

    @get_bond.register
    def _(self, atom_a_id: int, atom_b_id: int) -> Bond:
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
            Bond.from_dict(bond_data,topology.atoms)

        for angle_data in data["angles"]:
            Angle.from_dict(angle_data,topology.atoms)

        for dihedral_data in data["dihedrals"]:
            Dihedral.from_dict(dihedral_data,topology.atoms)

        for pair_data in data["pairs"]:
            Pair.from_dict(pair_data,topology.atoms)

        for exclusion_data in data["exclusions"]:
            Exclusion.from_dict(exclusion_data,topology.atoms)
        
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
        """
            Renumber the atoms.ids in the topology starting from a given starting number.
            This is used when extending a polymer with a new monomer.

            note: this is the only function that sets the formerly attribute of an atom
        """
        for atom in self.atoms:
            atom.formerly = atom.atom_id
            atom.atom_id = atom.atom_id + start 

    def add(self, topology):
        new_topology = copy.deepcopy(topology)
        self.atoms.extend(new_topology.atoms)

    def deduplicate(self):
        """
        Remove any duplicate bonds from the topology
        """
        # for each atom in the topology remove any duplicate bonds
        for atom in self.atoms:
            atom.deduplicate_bonds()

    def change_atom(self, old_atom: Atom, new_atom: Atom):
        """
        Change an atom in the topology to a new atom at the same position in the atom list
        """
        if not old_atom in self.atoms:
            raise ValueError("Atom not in topology")
        
        for bond in old_atom.bonds:
            bond.clone_bond_changing(old_atom, new_atom)
        
        # Find the index of old_atom
        index = self.atoms.index(old_atom)
        
        # Remove old_atom and insert new_atom at the same index
        self.atoms[index] = new_atom

        # remove pairs exclusions bonds (and associated angles, dihedrals) associated with the old atom
        old_atom.remove()

    def reverse(self) -> "Topology":
        copied_atoms = copy.deepcopy(self.atoms)
        reversed_atoms = copied_atoms[::-1]
        return Topology(reversed_atoms, self.preamble, self.molecule_type) 
    
    def contains_atom(self, candidate: Atom) -> bool:
        return any(atom for atom in self.atoms if atom == candidate)
    
    def contains_bond(self, candidate: Bond) -> bool:
        return any(bond for bond in self.bonds if bond == candidate)
    
    def __repr__(self) -> str:
        return f"Topology: ({len(self.atoms)} atoms, netcharge={self.netcharge})"
