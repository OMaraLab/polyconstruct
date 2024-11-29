from typing import List

class Junction:
    """Junctions are the polymerization sites of a monomer topology. 
    They are defined by a name (describing the PolymerJunction type) and a list of bonds.
    name: str - the name of the PolymerJunction type
    bonds: list - a list of bonds of that junction type
    """
    from polytop.atoms import Atom
    from polytop.bonds import Bond

    def __init__(self, remaining_atom : Atom, leaving_atom: Atom, name: str = None):
        if name is None:
            name = f"{remaining_atom.atom_name}-{leaving_atom.atom_name}"
        self.name = name
        self.remaining_atom = remaining_atom
        self.leaving_atom = leaving_atom
        if not leaving_atom in remaining_atom.bond_neighbours():
            raise ValueError("Junction location cannot be found")
        
    def named(self, newname : str) -> "Junction":
        self.name = newname
        return self
        
    def to_dict(self)->dict:
        return {
            "name": self.name,
            "remaining_atom": self.remaining_atom.atom_name,
            "leaving_atom": self.leaving_atom.atom_name
        }
    
    def second_remaining_atom(self)->Atom:
        return next((atom for atom in self.remaining_atom.bond_neighbours() if atom != self.leaving_atom), None) 

    def second_leaving_atom(self)->Atom:
        return next((atom for atom in self.leaving_atom.bond_neighbours() if atom != self.remaining_atom), None)

    @classmethod
    def from_dict(cls, data: dict, atoms: List[Atom]):
        name = data["name"]
        from polytop.atoms import Atom
        remaining_atom_name = data["remaining_atom"]
        remaining_atom = next(atom for atom in atoms if atom.atom_name == remaining_atom_name)
        leaving_atom_name = data["leaving_atom"]
        leaving_atom = next(atom for atom in atoms if atom.atom_name == leaving_atom_name)
        return cls(remaining_atom,leaving_atom,name)

    @classmethod
    def from_topology(cls, topology: "Topology", remaining_atom_name, leaving_atom_name, residue_id: int = None, name: str = None):
        remaining_atom = topology.get_atom(remaining_atom_name, residue_id)
        leaving_atom = topology.get_atom(leaving_atom_name, residue_id)
        return cls(remaining_atom, leaving_atom, name)

    def __repr__(self) -> str:
        return f"(\"{self.name}\":{self.remaining_atom.atom_name}-{self.leaving_atom.atom_name})"
    
class Junctions(list):
    def add(self, junction: Junction):
        self.append(junction)

    def get_junctions(self):
        return self

    def remove(self, junction: Junction):
        super().remove(junction)

    def named(self, name: str):
        return [junction for junction in self if junction.name == name]

    def to_dict(self):
        return [junction.to_dict() for junction in self]

    from polytop.atoms import Atom
    @classmethod
    def from_dict(cls, data: list, atoms: List[Atom]):
        junctions = cls()
        for junction_dict in data:
            junction = Junction.from_dict(junction_dict, atoms)
            junctions.add(junction)
        return junctions

    def __repr__(self) -> str:
        return f"({len(self)}) {','.join(j.__repr__() for j in self)}"