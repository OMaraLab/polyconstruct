from typing import List

class Junction:
    """Junctions are the polymerization sites of a monomer topology. 
    They are defined by a name (describing the PolymerJunction type) and a list of bonds.
    name: str - the name of the PolymerJunction type
    bonds: list - a list of bonds of that junction type
    """
    from .atoms import Atom
    from .bonds import Bond

    def __init__(self, monomer_atom : Atom, residue_atom: Atom, name: str = None):
        if name is None:
            name = f"{monomer_atom.atom_name}-{residue_atom.atom_name}"
        self.name = name
        self.monomer_atom = monomer_atom
        self.residue_atom = residue_atom
        if not residue_atom in monomer_atom.bond_neighbours():
            raise ValueError("Junction location cannot be found")
        
    def named(self, newname : str) -> "Junction":
        self.name = newname
        return self
        
    def to_dict(self)->dict:
        return {
            "name": self.name,
            "monomer_atom": self.monomer_atom.atom_name,
            "residue_atom": self.residue_atom.atom_name
        }

    @classmethod
    def from_dict(cls, data: dict, atoms: List[Atom]):
        name = data["name"]
        from polytop.atoms import Atom
        monomer_atom_name = data["monomer_atom"]
        monomer_atom = next(atom for atom in atoms if atom.atom_name == monomer_atom_name)
        residue_atom_name = data["residue_atom"]
        residue_atom = next(atom for atom in atoms if atom.atom_name == residue_atom_name)
        return cls(monomer_atom,residue_atom,name)

    @classmethod
    def from_topology(cls, topology: "Topology", monomer_atom_name, residue_atom_name, residue_id: int = None, name: str = None):
        monomer_atom = topology.get_atom(monomer_atom_name, residue_id)
        residue_atom = topology.get_atom(residue_atom_name, residue_id)
        return cls(monomer_atom, residue_atom, name)

    def __repr__(self) -> str:
        return f"(\"{self.name}\":{self.monomer_atom.atom_name}-{self.residue_atom.atom_name})"
    
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

    from .atoms import Atom
    @classmethod
    def from_dict(cls, data: list, atoms: List[Atom]):
        junctions = cls()
        for junction_dict in data:
            junction = Junction.from_dict(junction_dict, atoms)
            junctions.add(junction)
        return junctions

    # @classmethod
    # def from_list(cls, junctions_list: list["Junction"]):
    #     junctions = cls()
    #     for junction in junctions_list:
    #         junctions.add(junction)
    #     return junctions

    def __repr__(self) -> str:
        return f"({len(self)}) {','.join(j.__repr__() for j in self)}"
    
    
# class Junctions:
#     def __init__(self):
#         self.junctions = []
        
#     def __len__(self):
#         return len(self.junctions)
    
#     def __getitem__(self, index):
#         return self.junctions[index]

#     def add(self, junction: Junction):
#         self.junctions.append(junction)

#     def get_junctions(self):
#         return self.junctions
    
#     def remove(self, junction: Junction):
#         self.junctions.remove(junction)
    
#     def named(self, name: str):
#         return [junction for junction in self.junctions if junction.name == name]

#     def to_dict(self):
#         return [junction.to_dict() for junction in self.junctions]

#     @classmethod
#     def from_dict(cls, data: list, atoms: List[Atom]):
#         junctions = cls()
#         for junction_dict in data:
#             junction = Junction.from_dict(junction_dict, atoms)
#             junctions.add(junction)
#         return junctions
    
#     @classmethod
#     def from_list(cls, junctions_list: list["Junction"]):
#         junctions = cls()
#         for junction in junctions_list:
#             junctions.add(junction)
#         return junctions
#     def __repr__(self) -> str:
#         return f"Junctions({len(self.junctions)})"