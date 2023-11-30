from typing import List
from polytop.atoms import Atom
from polytop.bonds import Bond


class Junction:
    """Junctions are the polymerization sites of a monomer topology. 
    They are defined by a name (describing the PolymerJunction type) and a list of bonds.
    name: str - the name of the PolymerJunction type
    bonds: list - a list of bonds of that junction type
    """
    def __init__(self, name: str, location: Bond):
        self.name = name
        self.location = location
        if not location:
            raise ValueError("Junction location cannot be found")
        
    def to_dict(self)->dict:
        return {
            "name": self.name,
            "location": self.location.to_dict() if self.location else None
        }

    @classmethod
    def from_dict(cls, data: dict, atoms: List[Atom]):
        name = data["name"]
        location = Bond.from_dict(data["location"], atoms)
        return cls(name, location)
    

class Junctions:
    def __init__(self):
        self.junctions = []
        
    def __len__(self):
        return len(self.junctions)
    
    def __getitem__(self, index):
        return self.junctions[index]

    def add(self, junction: Junction):
        self.junctions.append(junction)

    def get_junctions(self):
        return self.junctions
    
    def remove(self, junction: Junction):
        self.junctions.remove(junction)
    
    def named(self, name: str):
        return [junction for junction in self.junctions if junction.name == name]

    def to_dict(self):
        return [junction.to_dict() for junction in self.junctions]

    @classmethod
    def from_dict(cls, data: list, atoms: List[Atom]):
        junctions = cls()
        for junction_dict in data:
            junction = Junction.from_dict(junction_dict, atoms)
            junctions.add(junction)
        return junctions
    
    @classmethod
    def from_list(cls, junctions_list: list["Junction"]):
        junctions = cls()
        for junction in junctions_list:
            junctions.add(junction)
        return junctions
    def __repr__(self) -> str:
        return f"Junctions({len(self.junctions)})"