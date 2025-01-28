from __future__ import annotations
import copy
from functools import singledispatchmethod
import re
import warnings
from typing import Dict, List, Optional, Tuple

from .Angles import Angle
from .Atoms import Atom
from .Bonds import Bond
from .Junction import Junction, Junctions
from .Dihedrals import Dihedral, Dihedral_type
from .Exclusions import Exclusion
from .Molecule_type import MoleculeType
from .Pairs import Pair


class Topology:
    """
    Represents the topology of a molecular system, including atoms, bonds,
    angles, and dihedrals.

    :param atoms: a list of atoms in the molecular system, defaults to None
    :type atoms: Optional[List[Atom]], optional
    :param preamble: a string containing the preamble of the topology file,
            defaults to None
    :type preamble: Optional[List[str]], optional
    :param molecule_type: the type of molecule represented by the topology,
            defaults to None
    :type molecule_type: Optional[MoleculeType], optional

    Attributes
    ----------
    atoms: List[Atom]
        A list of atoms in the molecular system.
    bonds: List[Bond]
        A list of bonds in the molecular system.
    angles: List[Angle]
        A list of angles in the molecular system.
    dihedrals: List[Dihedral]
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
        """
        Represents the topology of a molecular system, including atoms, bonds,
        angles, and dihedrals.

        :param atoms: a list of atoms in the molecular system, defaults to None
        :type atoms: Optional[List[Atom]], optional
        :param preamble: a string containing the preamble of the topology file,
                defaults to None
        :type preamble: Optional[List[str]], optional
        :param molecule_type: the type of molecule represented by the topology,
                defaults to None
        :type molecule_type: Optional[MoleculeType], optional
        """
        self.molecule_type = molecule_type
        self.atoms = atoms or []
        self.preamble = preamble or []
        self.title = "Unknown molecule"
        if self.preamble and self.preamble[1].startswith(';'):
            self.title = self.preamble[1].lstrip('; ')
        
    def copy(self) -> Topology:
        """
        A non-recursive alternative to deepcopy, which is used to create a new
        Topology with the same Atoms, Bonds, Angles, Dihedrals, etc. This
        function can copy large topology systems without the risk of stack
        overflow produced by a deepcopy algorithm.

        :return: a new copy of this Topology
        :rtype: Topology
        """
        new_atoms = []
        for atom in self.atoms:
            new_atoms.append(Atom.from_dict(atom.to_dict()))
        for pair in self.pairs:
            Pair.from_dict(pair.to_dict(),new_atoms)
        for exclusion in self.exclusions:
            Exclusion.from_dict(exclusion.to_dict(),new_atoms)
        for bond in self.bonds:
            Bond.from_dict(bond.to_dict(),new_atoms)
        for angle in self.angles:
            Angle.from_dict(angle.to_dict(),new_atoms)
        for dihedral in self.dihedrals:
            Dihedral.from_dict(dihedral.to_dict(),new_atoms)
        new_topology = Topology(
            atoms=new_atoms,
            preamble=self.preamble.copy(),
            molecule_type=copy.copy(self.molecule_type)
        )
        return new_topology
            
    @property
    def name(self) -> str:
        """
        Getter method to retrieve the Topology's name.

        :return: the name of the Topology, or "Unknown" if None.
        :rtype: str
        """
        return "Unknown" if self.molecule_type is None else self.molecule_type.name

    def __repr__(self) -> str:
        return f"Topology: {self.name} ({len(self.atoms)} atoms)"

    @property
    def pseudoatoms(self) -> dict[int: Atom]:
        """
        Getter method to retrieve the Topology's pseudo atoms.

        :return: a dictionary of the pseudo atoms in this Topology as values
                and their keys as their indexes.
        :rtype: dict[int: Atom]
        """
        return {atom.index: atom for atom in self.atoms if atom.is_virtual}

    @property
    def carbons(self) -> dict[int: Atom]:
        """
        Getter method to retrieve the Topology's carbon atoms.

        :return: a dictionary of the carbon atoms in this Topology as values
                and their keys as their indexes.
        :rtype: dict[int: Atom]
        """
        return {atom.index: atom for atom in self.atoms if atom.element == "C"}

    @property
    def bonds(self) -> List[Bond]:
        """
        Getter method to retrieve the Topology's Bonds.

        :return: a list of all of the Bonds in this Topology.
        :rtype: List[Bond]
        """
        bonds = []
        for atom in self.atoms:
            for bond in atom.bonds:
                if bond not in bonds:
                    bonds.append(bond)
        return bonds

    @property
    def angles(self) -> List[Angle]:
        """
        Getter method to retrieve the Topology's Angles.

        :return: a list of all of the Angles in this Topology.
        :rtype: List[Angle]
        """
        angles = []
        for atom in self.atoms:
            for bond in atom.bonds:
                for angle in bond.angles:
                    if angle not in angles:
                        angles.append(angle)
        return angles

    @property
    def dihedrals(self) -> List[Dihedral]:
        """
        Getter method to retrieve the Topology's Dihedrals.

        :return: a list of all of the Dihedrals in this Topology.
        :rtype: List[Dihedral]
        """
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
        """
        Getter method to retrieve the Topology's Pairs.

        :return: a list of all of the Pairs in this Topology.
        :rtype: List[Pair]
        """
        pairs = []
        for atom in self.atoms:
            for pair in atom.pairs:
                if pair not in pairs:
                    pairs.append(pair)
        return pairs

    @property
    def exclusions(self) -> List[Exclusion]:
        """
        Getter method to retrieve the Topology's Exclusions.

        :return: a list of all of the Exclusions in this Topology.
        :rtype: List[Exclusion]
        """
        exclusions = []
        for atom in self.atoms:
            for exclusion in atom.exclusions:
                if exclusion not in exclusions:
                    exclusions.append(exclusion)
        return exclusions
    
    @classmethod
    def numerically_order_oxygens(instance: Topology, section: str, line: str) -> str:
        """
        Preprocesses oxygen atoms with alphanumeric ordering to numeric 
        ordering. For example: OD becomes O4.

        :param instance: the Topology class being called.
        :type instance: Topology
        :param section: the itp file section, for example "molecule".
        :type section: str
        :param line: a line from the itp file.
        :type line: str
        :return: the processed line with oxygen atoms renamed numerically.
        :rtype: str
        """
        if section == 'atoms':
            new_line = re.sub(r'(\bO[A-Z]\b)', lambda match: 'O' + str(ord(match.group(1)[-1]) - ord('A') + 1), line) 
        else:
            new_line = line
        return new_line
    
    def _rearrange_atoms(self):
        """
        Rearrange/reorder atoms by their residue id and then by their rank.
        Additionally, ensure there are no 'gaps' in atom id values and all
        Atoms have a seperate charge group. Used before writing an ITP to
        ensure it is compliant with GROMACS requirements.
        """
        atoms = self.atoms.copy()
        
        atoms.sort(key=lambda atom: atom.residue_id)
        arranged_atoms = []
        for i in range(self.max_residue_id()+1):
            atom_subset = list(atom for atom in atoms if atom.residue_id==i)
            atom_subset.sort(key=lambda a: a.rank)
            arranged_atoms.extend(atom_subset)
        atoms = arranged_atoms

        for j, atom in enumerate(atoms):
            atom.atom_id = j+1
            atom.charge_group_num = j+1
        return atoms
    
    @classmethod
    def from_ITP(cls, file_path: str, preprocess=None, format: str = "gromos") -> Topology:
        """
        Class method to create a Topology object from an ITP file.

        Note: this function does not explicitly check that
        moleculetype.nrexcl = number of exclusions associated with atoms.

        :param file_path: the path to the GROMOS ITP file.
        :type file_path: str
        :param preprocess: function to preprocess the topology, that must
                accept the section and line as arguments in that order,
                defaults to None.
        :type preprocess: lambda function, optional
        :param format: The forcefield the ITP file is formatted as, options are
                "gromos", "amber", "opls" and "charmm"
        :type format: str, defaults to "gromos" for GROMOS forcefields.
        :return: the created Topology object.
        :rtype: Topology
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
                Angle.from_line(line, atoms)
            elif section == "pairs":
                Pair.from_line(line, atoms)
            elif section == "dihedrals":
                Dihedral.from_line(line, atoms, format=format)
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

    def to_ITP(self, file_path: str, format: str = "gromos"):
        """
        Write the Topology to a GROMACS ITP file of the desired forcefield format.

        :param file_path: the path to and the desired name of the GROMACS ITP
                file.
        :type file_path: str
        :param format: The forcefield the ITP file is formatted as, options are
                "gromos", "amber", "opls" and "charmm"
        :type format: str, defaults to "gromos" for GROMOS forcefields.
        """
        with open(file_path, "w") as f:
            if self.title:
                f.write(f";----------------------------TITLE -----------------------------------------------------------------------------------------\n")
                f.write(f"; {self.title}\n")
                f.write(f";\n")

            for line in self.preamble:
                if line.startswith(";"):
                    f.write(line + "\n")
                else:
                    f.write("; " + line + "\n")
            
            if self.molecule_type is not None:
                f.write("\n[ moleculetype ]\n")
                f.write(str(self.molecule_type) + "\n")

            f.write("[ atoms ]\n")
            self.atoms = self._rearrange_atoms()
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
                dihedral.dihedral_type.is_planar_constraint
                for dihedral in self.dihedrals
            ):
                f.write("\n[ dihedrals ]\n")
                f.write("; GROMOS improper dihedrals\n")
                for dihedral in self.dihedrals:
                    if dihedral.dihedral_type.is_planar_constraint:
                        f.write(str(dihedral) + "\n")

            f.write("\n[ dihedrals ]\n")
            for dihedral in self.dihedrals:
                if dihedral.dihedral_type.is_rotational_constraint:
                    f.write(str(dihedral) + "\n")

            f.write("\n[ exclusions ]\n")
            for exclusion in self.exclusions:
                f.write(str(exclusion) + "\n")
    @property 
    def residue_name(self) -> str:
        """
        Getter method to retrieve the Topology's residue name. The residue name
        of the entire Topology is set to that of the first Atom.

        :return: the residue name of the Topology.
        :rtype: str
        """
        return self.atoms[0].residue_name
    
    @residue_name.setter
    def residue_name(self, new_name):
        """
        Set the Topology's residue name to a desired string of 5 or less
        characters.

        :param new_name: the new residue name for the Topology.
        :type new_name: str
        :raises ValueError: if 'new_name' is longer than 5 characters.
        """
        if len(new_name) > 5:
            raise ValueError("Residue name must be 5 characters or less")
        for atom in self.atoms:
            atom.residue_name = new_name
        
    def atom_counts(self) -> Dict[str, int]:
        """
        Get the number of Atoms of each element in the Topology.

        :return: a dictionary of each element present as a key with the number
                of Atoms of this element as it's value.
        :rtype: Dict[str, int]
        """
        atom_counts = {}
        for atom in self.atoms:
            if atom.element not in atom_counts:
                atom_counts[atom.element] = 1
            else:
                atom_counts[atom.element] += 1
        return atom_counts
    
    def atom_elements (self) -> set[str]:
        """
        Retrieve a set of all atom elements in the Topology.

        :return: a set containing all of the elements present in this Topology.
        :rtype: set[str]
        """
        return set(atom.element for atom in self.atoms)
    
    def auto_rename_atoms(self):
        """
        Rename all atoms in the topology to their element symbol plus their
        index among atoms with the same element (ordered by atom.id).
        
        This compresses the namespace of atom names to the minimum possible,
        while preserving the element and order of atoms in the topology.

        For example, a carbon atom with the smallest id attribute out of
        all of the carbon atoms will be renamed to "C1".
        """
        for atom_symbol in self.atom_elements():
            self.reorder_atom_indexes(atom_symbol,1)
        
    def reorder_atom_indexes(self,atom_symbol: str,new_first_index: int):
        """
        Renumber the atom indexes of all atoms of a given symbol (i.e. element)
        to start from a given index. The original order of these atoms
        is maintained.

        :param atom_symbol: the element of atoms to renumber, e.g. "C" to
                renumber all of the carbon atoms
        :type atom_symbol: str
        :param new_first_index: the starting number for the selected atoms to
                be renumbered from, e.g. a value of '10' means that the first
                selected atom will now have an index of 10, and the next atom
                an index of 11 and so forth.
        :type new_first_index: int
        """
        atoms_with_this_symbol = [atom for atom in self.atoms if atom.element == atom_symbol]
        atoms_with_this_symbol.sort(key=lambda atom: atom.index)
        for atom in atoms_with_this_symbol:
            atom.index = new_first_index
            atom.element = atom_symbol
            new_first_index += 1

    def max_atom_index(self) -> dict[str, int]:
        """
        Get the maximum atom index for each atom symbol/element in the Topology. 

        :return: a dictionary with each element assigned a value representing
                the highest index of all atoms of that element.
        :rtype: dict[str, int]
        """
        max_atom_index = {}
        for atom in self.atoms:
            if atom.element not in max_atom_index:
                max_atom_index[atom.element] = atom.index
            else:
                max_atom_index[atom.element] = max(max_atom_index[atom.element], atom.index)
        return max_atom_index
            
    def junction(self, monomer_atom_name:str, residue_atom_name:str, name:str = None) -> Junction:
        """
        Additional method to create a Junction from atoms within this Topology.

        For example:
        top = Topology.from_ITP("file.itp")
        junction1 = top.junction("atomA", "atomB").named("name")

        Instead of:
        junction1 = Junction(top.get_atom("atomA"), top.get_atom("atomB"), name="name")

        :param monomer_atom_name: the name of the Atom which will remain with
                Topology after polymerisation, and will obtain a new bond.
        :type monomer_atom_name: str
        :param residue_atom_name: the name of the Atom which will be lost
                during polymerisation, analogous to the leaving atom in a
                chemical reaction.
        :type residue_atom_name: str
        :param name: the unique name of the Junction type, defaults to None
        :type name: str, optional
        :return: a new Junction
        :rtype: Junction
        """
        monomer_atom = self.get_atom(monomer_atom_name)
        residue_atom = self.get_atom(residue_atom_name)
        return Junction(monomer_atom, residue_atom, name)
    
    @property
    def residue_id(self) -> int:
        """
        Getter method to retrieve the Topology's residue id. The residue id
        of the entire Topology is set to that of the first Atom.

        :return: the residue id of the Topology.
        :rtype: int
        """
        return self.atoms[0].residue_id

    @residue_id.setter
    def residue_id(self, new_id: int):
        """
        Set the Topology's residue id to a desired value in the range of 1-99999.

        :param new_name: 
        :type new_name: str
        :raises ValueError: 

        :param new_id: the new residue id for the Topology.
        :type new_id: int
        :raises ValueError: if 'new_id' is less than 1.
        :raises ValueError: if 'new_id' is greater than 99999.
        """
        if new_id < 1:
            raise ValueError("Residue id must be greater than zero")
        if new_id > 99999:
            raise ValueError("Residue id must be less than 100000")
        for atom in self.atoms:
            atom.residue_id = new_id

    def max_residue_id(self) -> int:
        """
        Get the maximum residue id in the Topology. 

        :return: the highest residue id of all atoms in the Topology.
        :rtype: int
        """
        return max(atom.residue_id for atom in self.atoms)

    def renumber_residues(self, starting_from : int):
        """
        Renumber the residues by id in the topology starting from a given
        number. The order of the residues is maintained.

        :param starting_from: the starting number for the residues to
                be renumbered from, e.g. a value of '10' means that the first
                residue will now have an id of 10, and the next an id of 11
                and so forth.
        :type starting_from: int
        :raises ValueError: if 'starting_from' is less than 1.
        :raises ValueError: if any existing residue_id is less than 1.
        """
        if starting_from < 1:
            raise ValueError("Residue id to start from must be greater than zero")
        for atom in self.atoms:
            if atom.residue_id < 1:
                raise ValueError("All residue ids must be greater than zero")
        for atom in self.atoms:
            atom.residue_id += starting_from

    @property
    def netcharge(self) -> int:
        """
        Getter method to retrieve the Topology's net charge. The net charge is
        calculated as the sum of all of the Atom's partial charges.

        :return: the net charge of the Topology.
        :rtype: int
        """
        return sum(atom.partial_charge for atom in self.atoms)

    @netcharge.setter
    def netcharge(self, new_netcharge: float):
        """
        Adjusts the partial charges of all atoms by the same amount. To change
        the partial charge of each atom proportionally (i.e. to change the net
        charge without changing the charge distribution), use the function
        'proportional_charge_change' instead.

        :param new_netcharge: value to adjust the partial charges by.
        :type new_netcharge: float
        """
        old_netcharge = self.netcharge
        num_atoms_to_distribute = len(list(atom for atom in self.atoms if atom.residue_id == self.max_residue_id()))
        num_atoms_to_distribute += len(list(atom for atom in self.atoms if atom.residue_id == (self.max_residue_id()-1)))
        charge_difference_per_atom = (new_netcharge - old_netcharge) / num_atoms_to_distribute
        for atom in self.atoms:
            if atom.residue_id == self.max_residue_id() or atom.residue_id == (self.max_residue_id()-1):
                atom.partial_charge += charge_difference_per_atom

    def proportional_charge_change(self, new_netcharge: float):
        """
        Change the net charge of the molecular system by changing the partial 
        harge of each atom proportionally.

        :param new_netcharge: value to adjust the partial charges by.
        :type new_netcharge: float
        """
        old_netcharge = self.netcharge
        charge_difference = new_netcharge - old_netcharge
        total_absolute_charge = sum(abs(atom.partial_charge) for atom in self.atoms)
        
        for atom in self.atoms:
            atom.partial_charge += charge_difference * (abs(atom.partial_charge) / total_absolute_charge)

    def get_former_atom(self, former_atom_id: int) -> Atom:
        """
        Getter method to retrieve an Atom in the Topology from its 'formerly'
        attribute.

        :param former_atom_id: the former id (i.e. 'formerly' attribute) to retrieve an Atom.
        :type former_atom_id: int
        :return: an Atom from the Topology with a matching 'formerly' attribute.
        :rtype: Atom
        """
        return next((atom for atom in self.atoms if atom.formerly == former_atom_id), None)

    @singledispatchmethod
    def get_atom(self, atom_id: int) -> Atom:
        """
        Getter method to retrieve an Atom in the Topology from its id.

        :param atom_id: the id of the desired Atom.
        :type atom_id: int
        :return: the Atom corresponding to the provided 'atom_id'.
        :rtype: Atom
        """
        atom = next((atom for atom in self.atoms if atom.atom_id == atom_id), None)
        return atom

    
    @get_atom.register
    def _(self, atom_name: str, residue_id : int = None) -> Atom:
        """
        Getter method to retrieve an Atom in the Topology from its name, and a
        residue_id if the name is not unique.

        :param atom_name: the name of the desired Atom.
        :type atom_name: str
        :param residue_id: the residue_id of the desired Atom, defaults to None.
        :type residue_id: int, optional
        :raises ValueError: if there is no Atom in the Topology with this name.
        :raises ValueError: if there is more than one Atom in the Topology with
                this name.
        :raises ValueError: if there is no Atom in the Topology with this name
                and residue_id.
        :return: the corresponding Atom.
        :rtype: Atom
        """
        atoms = [atom for atom in self.atoms if atom.atom_name == atom_name]
        if len(atoms) == 0:
            raise ValueError(f"No atom with name {atom_name} in the topology.")
        if residue_id is None:
            if len(atoms) > 1:
                raise ValueError(f"Multiple atoms with id {atom_name} in the topology. Please specify the residue id.")
            return atoms[0]
        else:
            atom = next((atom for atom in atoms if atom.residue_id == residue_id), None)
            if atom is None:
                raise ValueError(f"No atom with name {atom_name} and residue id {residue_id} in the topology.")
            return atom

    @singledispatchmethod
    def get_bond(self, atom_a: Atom, atom_b: Atom) -> Bond:
        """
        Getter method to retrieve a Bond in the Topology from the two Atoms
        that share it.

        :param atom_a: The first Atom involved in the bond.
        :type atom_a: Atom
        :param atom_b: The second Atom involved in the bond.
        :type atom_b: Atom
        :return: the corresponding Bond.
        :rtype: Bond
        """
        return Bond.from_atoms(atom_a, atom_b)

    @get_bond.register
    def _(self, atom_a_id: int, atom_b_id: int) -> Bond:
        """
        Getter method to retrieve a Bond in the Topology from the id's of the
        two Atoms that share it.

        :param atom_a_id: The id of the first Atom involved in the bond.
        :type atom_a_id: int
        :param atom_b_id: The id of the second Atom involved in the bond.
        :type atom_b_id: int
        :return: the corresponding Bond.
        :rtype: Bond
        """
        atom_a = self.get_atom(atom_a_id)
        atom_b = self.get_atom(atom_b_id)
        return Bond.from_atoms(atom_a, atom_b)

    @get_bond.register
    def _(self, atom_a_name: str, atom_b_name: str) -> Bond:
        """
        Getter method to retrieve a Bond in the Topology from the names of the
        two Atoms that share it.

        :param atom_a_name: The name of the first Atom involved in the bond.
        :type atom_a_name: str
        :param atom_b_name: The name of the second Atom involved in the bond.
        :type atom_b_name: str
        :return: the corresponding Bond.
        :rtype: Bond
        """
        atom_a = self.get_atom(atom_a_name)
        atom_b = self.get_atom(atom_b_name)
        return Bond.from_atoms(atom_a, atom_b)

    @singledispatchmethod
    def get_pair(self, atom_a: Atom, atom_b: Atom) -> Pair:
        """
        Getter method to retrieve a Pair in the Topology from the
        two Atoms that share it.

        :param atom_a: The first Atom involved in the pair.
        :type atom_a: Atom
        :param atom_b: The second Atom involved in the pair.
        :type atom_b: Atom
        :return: the corresponding Pair.
        :rtype: Pair
        """
        return Pair.from_atoms(atom_a, atom_b)

    @get_pair.register
    def _(self, atom_a_id: int, atom_b_id: int) -> Pair:
        """
        Getter method to retrieve a Pair in the Topology from the id's of the
        two Atoms that share it.

        :param atom_a_id: The id of the first Atom involved in the pair.
        :type atom_a_id: int
        :param atom_b_id: The id of the second Atom involved in the pair.
        :type atom_b_id: int
        :return: the corresponding Pair.
        :rtype: Pair
        """
        atom_a = self.get_atom(atom_a_id)
        atom_b = self.get_atom(atom_b_id)
        return Pair.from_atoms(atom_a, atom_b)

    @get_pair.register
    def _(self, atom_a_name: str, atom_b_name: str) -> Pair:
        """
        Getter method to retrieve a Pair in the Topology from the names of the
        two Atoms that share it.

        :param atom_a_name: The name of the first Atom involved in the pair.
        :type atom_a_name: str
        :param atom_b_name: The name of the second Atom involved in the pair.
        :type atom_b_name: str
        :return: the corresponding Pair.
        :rtype: Pair
        """
        atom_a = self.get_atom(atom_a_name)
        atom_b = self.get_atom(atom_b_name)
        return Pair.from_atoms(atom_a, atom_b)

    @singledispatchmethod
    def get_angle(self, atom_a_id: int, atom_b_id: int, atom_c_id: int) -> Angle:
        """
        Getter method to retrieve an Angle in the Topology from the id's of the
        three Atoms that share it.

        :param atom_a_id: The id of the first Atom involved in the angle.
        :type atom_a_id: int
        :param atom_b_id: The id of the second Atom involved in the angle.
        :type atom_b_id: int
        :param atom_c_id: The id of the third Atom involved in the angle.
        :type atom_c_id: int
        :return: the corresponding Angle.
        :rtype: Angle
        """
        atom_a = self.get_atom(atom_a_id)
        atom_b = self.get_atom(atom_b_id)
        atom_c = self.get_atom(atom_c_id)
        return Angle.from_atoms(atom_a, atom_b, atom_c)

    @get_angle.register
    def _(self, atom_a_name: str, atom_b_name: str, atom_c_name: str) -> Angle:
        """
        Getter method to retrieve an Angle in the Topology from the names of the
        three Atoms that share it.

        :param atom_a_name: The name of the first Atom involved in the angle.
        :type atom_a_name: str
        :param atom_b_name: The name of the second Atom involved in the angle.
        :type atom_b_name: str
        :param atom_c_name: The name of the third Atom involved in the angle.
        :type atom_c_name: str
        :return: the corresponding Angle.
        :rtype: Angle
        """
        atom_a = self.get_atom(atom_a_name)
        atom_b = self.get_atom(atom_b_name)
        atom_c = self.get_atom(atom_c_name)
        return Angle.from_atoms(atom_a, atom_b, atom_c)

    @singledispatchmethod
    def get_dihedral(
        self, atom_a_id: int, atom_b_id: int, atom_c_id: int, atom_d_id: int
    ) -> Dihedral:
        """
        Getter method to retrieve a Dihedral in the Topology from the id's of
        the four Atoms that share it.

        :param atom_a_id: The id of the first Atom involved in the dihedral.
        :type atom_a_id: int
        :param atom_b_id: The id of the second Atom involved in the dihedral.
        :type atom_b_id: int
        :param atom_c_id: The id of the third Atom involved in the dihedral.
        :type atom_c_id: int
        :param atom_d_id: The id of the fourth Atom involved in the dihedral.
        :type atom_d_id: int
        :return: the corresponding Dihedral.
        :rtype: Dihedral
        """
        atom_a = self.get_atom(atom_a_id)
        atom_b = self.get_atom(atom_b_id)
        atom_c = self.get_atom(atom_c_id)
        atom_d = self.get_atom(atom_d_id)
        return Dihedral.from_atoms(atom_a, atom_b, atom_c, atom_d)

    @get_dihedral.register
    def _(self, atom_a_name: str, atom_b_name: str, atom_c_name: str, atom_d_name: str) -> Dihedral:
        """
        Getter method to retrieve a Dihedral in the Topology from the names of
        the four Atoms that share it.

        :param atom_a_name: The name of the first Atom involved in the dihedral.
        :type atom_a_name: str
        :param atom_b_name: The name of the second Atom involved in the dihedral.
        :type atom_b_name: str
        :param atom_c_name: The name of the third Atom involved in the dihedral.
        :type atom_c_name: str
        :param atom_d_name: The name of the fourth Atom involved in the dihedral.
        :type atom_d_name: str
        :return: the corresponding Dihedral.
        :rtype: Dihedral
        """
        atom_a = self.get_atom(atom_a_name)
        atom_b = self.get_atom(atom_b_name)
        atom_c = self.get_atom(atom_c_name)
        atom_d = self.get_atom(atom_d_name)
        return Dihedral.from_atoms(atom_a, atom_b, atom_c, atom_d)

    def add_atom(self, atom: Atom):
        """
        Add a new Atom to the Topology. The Atom will be appended to the end of
        the list of Atoms.

        :param atom: a new Atom to add to the Topology.
        :type atom: Atom
        """
        self.atoms.append(atom)

    def to_dict(self) -> dict:
        """
        Convert this Topology to a dictionary representation.

        The structure of the dictionary is as below:
        {"atoms": [atom.to_dict() for atom in self.atoms],
        "bonds": [bond.to_dict() for bond in self.bonds],
        "angles": [angle.to_dict() for angle in self.angles],
        "dihedrals": [dihedral.to_dict() for dihedral in self.dihedrals],
        "pairs": [pair.to_dict() for pair in self.pairs],
        "exclusions": [exclusion.to_dict() for exclusion in self.exclusions],
        "preamble": self.preamble,
        "moleculetype": self.molecule_type.to_dict()}

        :return: a dictionary representation of this Topology's attributes.
        :rtype: dict
        """
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
    def from_dict(cls, data: dict) -> Topology:
        """
        Create a new Topology from a dictionary, such as that created with
        Topology.to_dict().

        The structure of the dictionary is as below:
        {"atoms": [atom.to_dict() for atom in self.atoms],
        "bonds": [bond.to_dict() for bond in self.bonds],
        "angles": [angle.to_dict() for angle in self.angles],
        "dihedrals": [dihedral.to_dict() for dihedral in self.dihedrals],
        "pairs": [pair.to_dict() for pair in self.pairs],
        "exclusions": [exclusion.to_dict() for exclusion in self.exclusions],
        "preamble": self.preamble,
        "moleculetype": self.molecule_type.to_dict()}

        :param data: dictionary containing data to make a Topology, generate
                with 'to_dict()'.
        :type data: dict
        :return: a new Topology
        :rtype: Topology
        """
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
        """
        Remove an Atom and all Bonds, Angles and Dihedrals associated with it.
        If the atom is not present, exception errors are ignored.

        :param atom: the Atom to be removed from the Topology
        :type atom: Atom
        """
        try:
            atom.remove()
            self.atoms.remove(atom)        
        except ValueError:
            pass

    def reorder_atoms(self):
        """
        Renumber all of the atoms atom_id attributes starting from 1. The
        existing order is maintained.
        """
        for i in range(len(self.atoms)):
            self.atoms[i].atom_id = i + 1  # note atom id's are 1 indexed
            
    def renumber_atoms(self, start: int):
        """
        Renumber all of the atoms atom_id attributes to their current value
        plus the value of 'start'. The existing order is maintained.
        
        This is used when extending a polymer with a new monomer.

        Note: this is the only function that sets the 'formerly' attribute of
            an atom

        :param start: the value to offset the Atom's 'atom_id' attributes by.
        :type start: int
        """
        for atom in self.atoms:
            atom.formerly = atom.atom_id
            atom.atom_id = atom.atom_id + start 

    def set_former_ids(self):
        """
        Sets the 'formerly' attribute of all of the Atoms in the Topology to
        their current atom_id.
        """
        for atom in self.atoms:
            atom.formerly = atom.atom_id

    def clear_former_ids(self):
        """
        Sets the 'formerly' attribute of all of the Atoms in the Topology
        to None.
        """
        for atom in self.atoms:
            atom.formerly = None
            
    def __add__(self, other: Topology) -> Topology:
        this_topology = self.copy()
        this_topology.add(other)
        return this_topology
    
    def add(self, topology: Topology):
        """
        Add a new Topology 'topology' to this Topology.

        :param topology: a Topology to be added to this Topology.
        :type topology: Topology
        """
        new_topology = topology.copy()
        # set the residue id of the new atoms to the maximum residue id
        # in the current topology plus 1
        new_residue_id = self.max_residue_id() + 1
        for atom in new_topology.atoms:
            atom.residue_id = new_residue_id
        self.atoms.extend(new_topology.atoms)
        self.reorder_atoms()  # reorder atoms so the atom_ids are correct

    def deduplicate(self):
        """
        Remove any duplicate bonds from the topology.
        """
        # for each atom in the topology remove any duplicate bonds
        for atom in self.atoms:
            atom.deduplicate_bonds()

    def change_atom(self, old_atom: Atom, new_atom: Atom):
        """
        Change an Atom in the Topology to a new Atom at the same position
        in the atom list. Bonds, Angles, Dihedrals, etc. are also updated.

        :param old_atom: the Atom to swap out of the Atom list.
        :type old_atom: Atom
        :param new_atom: the Atom to swap into the Atom list.
        :type new_atom: Atom
        :raises ValueError: if 'old_atom' is not in the Topology.
        """
        if not old_atom in self.atoms:
            raise ValueError("Atom not in topology")

        # change all the bonds
        for bond in old_atom.bonds:
            bond.clone_bond_changing(old_atom, new_atom)
        
        # change all the angles
        for bond in old_atom.bonds:
            for angle in bond.angles:
                angle.clone_angle_changing(old_atom, new_atom)

        # change all the dihedrals
        for bond in old_atom.bonds:
            for angle in bond.angles:
                for dihedral in angle.dihedrals:
                    dihedral.clone_dihedral_changing(old_atom, new_atom)

        # Find the index of old_atom
        index = self.atoms.index(old_atom)

        self.atoms[index] = new_atom
        
        # remove pairs exclusions bonds (and associated angles, dihedrals) associated with the old atom
        old_atom.remove()

    def reverse(self) -> Topology:
        """
        Create a new Topology that is a copy of this Topology, except that the
        atoms are in the opposite order.

        :return: the new, reverse order Topology.
        :rtype: Topology
        """
        copied_atoms = self.copy().atoms
        reversed_atoms = copied_atoms[::-1]
        return Topology(reversed_atoms, self.preamble, self.molecule_type) 
    
    def contains_atom(self, candidate: Atom) -> bool:
        """
        Check if the Topology contains the given Atom.

        :param candidate: the Atom that is being checked to see if it is
                in the Topology.
        :type candidate: Atom
        :return: True if the Atom is in the Topology, otherwise False
        :rtype: bool
        """
        return any(atom for atom in self.atoms if atom == candidate)
    
    def contains_bond(self, candidate: Bond) -> bool:
        """
        Check if the Topology contains the given Bond.

        :param candidate: the Bond that is being checked to see if it is
                in the Topology.
        :type candidate: Bond
        :return: True if the Bond is in the Topology, otherwise False
        :rtype: bool
        """
        return any(bond for bond in self.bonds if bond == candidate)
    
    def __repr__(self) -> str:
        atom_names = [atom.atom_name for atom in self.atoms]
        atom_names_str = ",".join(atom_names)
        return f"({len(self.atoms)}) [{atom_names_str}] netcharge={self.netcharge}"
