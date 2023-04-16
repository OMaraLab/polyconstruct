import json
from polytop.molecule_type import MoleculeType

def test_serialization():
    moltype = MoleculeType("ARG", 3)  
    with open("tests/samples/moleculetype.json", "w") as f:
        json.dump(moltype.to_dict(), f)
    with open("tests/samples/moleculetype.json", "r") as f:
        new_moleculetype = MoleculeType.from_dict(json.load(f))
    assert moltype.name == new_moleculetype.name
    assert moltype.nrexcl == new_moleculetype.nrexcl
    
    