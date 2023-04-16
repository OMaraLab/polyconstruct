class MoleculeType:
    def __init__(self, name: str, nrexcl: int):
        self.name = name
        self.nrexcl = nrexcl

    @classmethod
    def from_line(cls, line: str) -> "MoleculeType":
        parts = line.split()
        name = parts[0]
        nrexcl = int(parts[1])
        return cls(name, nrexcl)

    def __str__(self) -> str:
        return f"{self.name:5} {self.nrexcl:5}"
    
    def __repr__(self) -> str:
        return f"{self.name:5} {self.nrexcl:5}"
    
    def __eq__(self, other: "MoleculeType") -> bool:
        return self.name == other.name and self.nrexcl == other.nrexcl

    def to_dict(self):
        moltype_dict = {
            "name": self.name,
            "nrexcl": self.nrexcl,
        }
        return moltype_dict

    @classmethod
    def from_dict(cls, data):
        return cls(
            name=data["name"],
            nrexcl=data["nrexcl"],
        )
        