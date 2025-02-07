from __future__ import annotations

from typing import Dict, List, Union

class Bond:
    """
    Represents a bond between two atoms in a molecular system.

    :param atom_a: The first atom involved in the bond.
    :type atom_a: Atom
    :param atom_b: The second atom involved in the bond.
    :type atom_b: Atom
    :param bond_type: The type of the bond (e.g., single, double, triple).
    :type bond_type: int
    :param bond_length: The length of the bond.
    :type bond_length: float
    :param force_constant: The force constant associated with the bond.
    :type force_constant: float
    :param order: The bond order, defaults to 1 (single bond)
    :type order: int, optional
    :raises ValueError: if either atom_a or atom_b are None
    """
    def __init__(
        self,
        atom_a: "Atom",
        atom_b: "Atom",
        bond_type: int,
        bond_length: float,
        force_constant: float,
        order: int = 1,
        format: str = "gromos"
    ) -> None:
        """
        Represents a bond between two atoms in a molecular system.

        :param atom_a: The first atom involved in the bond.
        :type atom_a: Atom
        :param atom_b: The second atom involved in the bond.
        :type atom_b: Atom
        :param bond_type: The type of the bond (e.g., single, double, triple).
        :type bond_type: int
        :param bond_length: The length of the bond.
        :type bond_length: float
        :param force_constant: The force constant associated with the bond.
        :type force_constant: float
        :param order: The bond order, defaults to 1 (single bond)
        :type order: int, optional
        :param format: The forcefield the ITP file is formatted as, options are
                "gromos", "amber", "opls" and "charmm"
        :type format: str, defaults to "gromos" for GROMOS forcefields.
        :raises ValueError: if either atom_a or atom_b are None
        """
        if atom_a is None or atom_b is None:
            raise ValueError("Bond must have two atoms")
        self.atom_a = atom_a
        self.atom_b = atom_b
        if self.atom_a.atom_id > self.atom_b.atom_id:  # keep atom order ascending
            self.atom_a, self.atom_b = self.atom_b, self.atom_a
        self.bond_type = bond_type
        self.bond_length = bond_length
        self.force_constant = force_constant
        atom_a.bonds.add(self)
        atom_b.bonds.add(self)
        self.order = order
        self.angles = set()
        self.format = format

    @classmethod
    def from_line(cls, line: str, atoms: List["Atom"], format: str = "gromos") -> Bond:
        """
        Class method to construct Bond from the line of an ITP file and a list
        of all Atom's present in the topology.

        :param line: the ITP file line
        :type line: str
        :param atoms: list of all Atoms in the Topology 
        :type atoms: List[Atom]
        :param format: The forcefield the ITP file is formatted as, options are
                "gromos", "amber", "opls" and "charmm"
        :type format: str, defaults to "gromos" for GROMOS forcefields.
        :return: the new Bond
        :rtype: Bond
        """
        parts = line.split()
        atom_a = atoms[int(parts[0]) - 1]
        atom_b = atoms[int(parts[1]) - 1]
        bond_type = int(parts[2])
        if format=="charmm":
            bond_length = 0 # 0 instead of None enables extend to work, is not saved to output
            force_constant = 0
        else:
            bond_length = float(parts[3])
            force_constant = float(parts[4])

        return cls(atom_a, atom_b, bond_type, bond_length, force_constant, format=format)

    @staticmethod
    def from_atoms(atom_a: "Atom", atom_b: "Atom", find_empty: bool = False) -> Bond:
        """
        Class method to find and return Bond from between two Atoms. 

        :param atom_a: The first atom involved in the angle.
        :type atom_a: Atom
        :param atom_b: The second atom in the angle.
        :type atom_b: Atom
        :param find_empty: Optional argument used when de-duplicating bonds to
                ensure the bond without angles associated is returned to delete.
        :type find_empty: bool, defaults to False
        :return: a Bond between these Atoms, or None if either atom_a or atom_b
                are None or there is not a bond between them.
        :rtype: Bond
        """
        if atom_a is None or atom_b is None:
            return None

        send = list(bond for bond in atom_a.bonds if bond.atom_b == atom_b or bond.atom_a == atom_b)
        # need to prevent ambiguity in selecting bonds, so that only one duplicate
        # bond created during polymer extension gets all of the Angles and Dihedrals,
        # and the other is removed during bond de-duplication
        if len(send) > 1 and find_empty==False:
            # duplicate bond, return the bond that has angles already to give it more
            has_angles = list(bond for bond in send if len(bond.angles)>=1 and bond.angles!=None and bond.angles!=set())[0]
            return has_angles
        elif len(send) > 1 and find_empty==True:
            # duplicate bond, return the bond with no angles to delete
            no_angles = list(bond for bond in send if len(bond.angles)==0 or bond.angles==None or bond.angles==set())[0]
            return no_angles
        elif len(send)==1:
            # only 1 bond for this atom pair, return it
            return send[0]
        else:
            # this atom pair does not have a shared bond
            return None

    def contains_atom(self, atom: "Atom") -> bool:
        """
        Check if this Bond contains a given atom.

        :param atom: the Atom you wish to check if it is in this Bond or not
        :type atom: Atom
        :return: True if the Bond contains the given "Atom", or False if not.
        :rtype: bool
        """
        return atom in [self.atom_a, self.atom_b]

    def clone_bond_changing(self, from_atom: "Atom", to_atom: "Atom") -> Bond:
        """
        Clone the bond, changing the atom that is being replaced. Used during
        the polymer.extend() algorithm to copy and modify bonds where a new
        Monomer is joined to the Polymer.

        :param from_atom: the outgoing "Atom", to be replaced
        :type from_atom: Atom
        :param to_atom: the incoming "Atom", will replace the position of the
                outgoing Atom in this Bond
        :type to_atom: Atom
        :raises ValueError: if 'from_atom' is not in the Bond
        :return: the new, modified Bond
        :rtype: Bond
        """
        if self.atom_a == from_atom: # first atom is being replaced
            new_bond = Bond(to_atom, self.atom_b, self.bond_type, self.bond_length, self.force_constant, self.order, self.format)
        elif self.atom_b == from_atom: # second atom is being replaced
            new_bond = Bond(self.atom_a, to_atom, self.bond_type, self.bond_length, self.force_constant, self.order, self.format)
        else:
            raise ValueError(f"Atom {from_atom} is not in bond {self}")
        return new_bond

    def other_atom(self, atom: "Atom")-> "Atom":
        """
        Check if the given Atom is in this Angle and return a list of the other
        atoms present in this Angle (i.e. discluding 'atom').

        :param atom: the Atom you wish to check if it is in this Bond, and if
                so, which other atom it is bonded to
        :type atom: Atom
        :raises ValueError: if 'atom' is not in this Bond
        :return: the other Atom in this Bond
        :rtype: Atom
        """
        if atom == self.atom_a:
            return self.atom_b
        elif atom == self.atom_b:
            return self.atom_a
        else:
            raise ValueError(f"Atom {atom} is not in bond {self}")
    
    def LHS(self) -> set["Atom"]:
        """
        List of all atoms in the left-hand side of the bond.

        :return: set of Atoms located on the left-hand side of this bond
        :rtype: set[Atom]
        """
        LHS_atoms = set()
        def traverse(atom: "Atom"):
            if atom != self.atom_b:
                LHS_atoms.add(atom)
                neighbours = atom.bond_neighbours()
                for neighbour in list(neighbours):
                    if neighbour not in LHS_atoms and neighbour != self.atom_b:
                        traverse(neighbour)

        traverse(self.atom_a)
        return LHS_atoms
    
    def RHS(self) -> set["Atom"]:
        """
        List of all atoms in the right-hand side of the bond.

        :return: set of Atoms located on the right-hand side of this bond
        :rtype: set[Atom]
        """
        RHS_atoms = set()
        def traverse(atom: "Atom"):
            if atom != self.atom_a:
                RHS_atoms.add(atom)
                neighbours = atom.bond_neighbours()
                for neighbour in list(neighbours):
                    if neighbour not in RHS_atoms and neighbour != self.atom_a:
                        traverse(neighbour)
        traverse(self.atom_b)
        return RHS_atoms
        
    def remove(self):
        """
        Delete self from all related Angles and the two Atoms. Used to cleanup
        and remove attributes during Polymer.extend().
        """
        while self.angles:
            self.angles.pop().remove()
        if self in self.atom_a.bonds:
            self.atom_a.bonds.remove(self)
        if self in self.atom_b.bonds:
            self.atom_b.bonds.remove(self)
                
    def __str__(self):
        if self.format == "charmm":
            return f"{self.atom_a.atom_id:>5} {self.atom_b.atom_id:>5} {self.bond_type:>5}"
        else:
            return f"{self.atom_a.atom_id:>5} {self.atom_b.atom_id:>5} {self.bond_type:>5} {self.bond_length:>10.4f} {self.force_constant:.4e}"

    def __repr__(self) -> str:
        if self.order == 1:
            return f"Bond({self.atom_a.atom_id} - {self.atom_b.atom_id})"
        elif self.order == 2:
            return f"Bond({self.atom_a.atom_id} = {self.atom_b.atom_id})"
        elif self.order == 3:
            return f"Bond({self.atom_a.atom_id} ≡ {self.atom_b.atom_id})"
        elif self.order == 0:
            return f"Bond({self.atom_a.atom_id} | {self.atom_b.atom_id})"
        else:
            return f"Bond({self.atom_a.atom_id} {self.atom_b.atom_id})"

    def to_dict(self) -> dict:
        """
        Convert this Bond to a dictionary representation.

        The structure of the dictionary is as below:
        {"atom_a": self.atom_a.atom_id,
        "atom_b": self.atom_b.atom_id,
        "bond_type": self.bond_type,
        "bond_length": self.bond_length,
        "force_constant": self.force_constant,
        "order": self.order}

        :return: a dictionary containing the id's of its Atoms and other
                attributes of this Bond.
        :rtype: dict
        """
        return {
            "atom_a": self.atom_a.atom_id,
            "atom_b": self.atom_b.atom_id,
            "bond_type": self.bond_type,
            "bond_length": self.bond_length,
            "force_constant": self.force_constant,
            "order": self.order,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Union[int, float]], atoms: List["Atom"]) -> Bond:
        """
        Create a new Bond from a dictionary, such as that created with
        Bond.to_dict(). Will retrieve an existing Bond if
        it already exists between these Atoms.

        The structure of the dictionary is as below:
        {"atom_a": self.atom_a.atom_id,
        "atom_b": self.atom_b.atom_id,
        "bond_type": self.bond_type,
        "bond_length": self.bond_length,
        "force_constant": self.force_constant,
        "order": self.order}

        :param data: dictionary containing data to make a Bond, generate
                with 'to_dict()'.
        :type data: Dict[str, Union[int, float]]
        :param atoms: list of Atoms. The list may contain more than 3 atoms, as
                long as the id's of the three atoms specified in the data dict
                are present.
        :type atoms: List[Atom]
        :return: a new Bond
        :rtype: Bond
        """
        atom_a = next((atom for atom in atoms if atom.atom_id == data['atom_a']),None)
        atom_b = next((atom for atom in atoms if atom.atom_id == data['atom_b']), None)
        # check for existing bond
        existing_bond = Bond.from_atoms(atom_a,atom_b)
        if existing_bond:
            return existing_bond
        else:
            return cls(
                atom_a = atom_a,
                atom_b = atom_b,
                bond_type=data["bond_type"],
                bond_length=data["bond_length"],
                force_constant=data["force_constant"],
                order=data["order"],
            )
        
    # def __eq__(self, __value: object) -> bool:
    #     if isinstance(__value, Bond):
    #         return self.atom_a.atom_id == __value.atom_a.atom_id and self.atom_b.atom_id == __value.atom_b.atom_id
    #     else:
    #         return False
        
    # def __hash__(self) -> int:
    #     hash_value = hash((self.atom_a, self.atom_b))
    #     return hash_value