from __future__ import annotations

from typing import Dict, List, Union

class Atom:
    ...

class Bond:
    """
    Represents a bond between two atoms in a molecular system.

    Attributes
    ----------
    atom_a : Atom
        The first atom involved in the bond.
    atom_b : Atom
        The second atom involved in the bond.
    bond_type : int
        The type of the bond (e.g., single, double, triple).
    bond_length : float
        The length of the bond.
    force_constant : float
        The force constant associated with the bond.
    bond_order : int, optional
        The bond order, default is 1 (single bond).

    """
    def __init__(
        self,
        atom_a: Atom,
        atom_b: Atom,
        bond_type: int,
        bond_length: float,
        force_constant: float,
        order: int = 1,
    ) -> None:
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

    @classmethod
    def from_line(cls, line: str, atoms: List["Atom"]):
        parts = line.split()
        atom_a = atoms[int(parts[0]) - 1]
        atom_b = atoms[int(parts[1]) - 1]
        bond_type = int(parts[2])
        bond_length = float(parts[3])
        force_constant = float(parts[4])

        return cls(atom_a, atom_b, bond_type, bond_length, force_constant)

    @staticmethod
    def from_atoms(atom_a: Atom, atom_b: Atom):
        if atom_a is None or atom_b is None:
            return None
        return next(
            (
                bond
                for bond in atom_a.bonds
                if bond.atom_b == atom_b or bond.atom_a == atom_b
            ),
            None,
        )

    def contains_atom(self, atom: Atom) -> bool:
        return atom in [self.atom_a, self.atom_b]

    def clone_bond_changing(self, from_atom: Atom, to_atom: Atom):
        """ Clone the bond, changing the atom that is being replaced """
        if self.atom_a == from_atom: # first atom is being replaced
            new_bond = Bond(to_atom, self.atom_b, self.bond_type, self.bond_length, self.force_constant, self.order)
        elif self.atom_b == from_atom: # second atom is being replaced
            new_bond = Bond(self.atom_a, to_atom, self.bond_type, self.bond_length, self.force_constant, self.order)
        else:
            raise ValueError(f"Atom {from_atom} is not in bond {self}")
        return new_bond

    def references_atom(self, atom: Atom) -> bool:
        return atom in [self.atom_a, self.atom_b]

    def other_atom(self, atom: Atom)-> Atom:
        if atom == self.atom_a:
            return self.atom_b
        elif atom == self.atom_b:
            return self.atom_a
        else:
            raise ValueError(f"Atom {atom} is not in bond {self}")
    
    def LHS(self) -> set['Atom']:
        # list all atoms in the LHS of the bond
        LHS_atoms = set()
        def traverse(atom: Atom):
            if atom != self.atom_b:
                LHS_atoms.add(atom)
                neighbours = atom.bond_neighbours()
                for neighbour in list(neighbours):
                    if neighbour not in LHS_atoms and neighbour != self.atom_b:
                        traverse(neighbour)

        traverse(self.atom_a)
        return LHS_atoms
    
    def RHS(self) -> set['Atom']:
        # list all atoms in the RHS of the bond
        RHS_atoms = set()
        def traverse(atom: Atom):
            if atom != self.atom_a:
                RHS_atoms.add(atom)
                neighbours = atom.bond_neighbours()
                for neighbour in list(neighbours):
                    if neighbour not in RHS_atoms and neighbour != self.atom_a:
                        traverse(neighbour)
        traverse(self.atom_b)
        return RHS_atoms
        
    def remove(self):
        while self.angles:
            self.angles.pop().remove()
        if self in self.atom_a.bonds:
            self.atom_a.bonds.remove(self)
        if self in self.atom_b.bonds:
            self.atom_b.bonds.remove(self)
                
    def __str__(self):
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

    def to_dict(self):
        return {
            "atom_a": self.atom_a.atom_id,
            "atom_b": self.atom_b.atom_id,
            "bond_type": self.bond_type,
            "bond_length": self.bond_length,
            "force_constant": self.force_constant,
            "order": self.order,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Union[int, float]], atoms: List[Atom]):
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