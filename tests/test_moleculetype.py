import json
from pathlib import Path
from polytop.Molecule_type import MoleculeType

def test_serialization(output_dir: Path):
    moltype = MoleculeType("ARG", 3) 
    
    # Write to the file
    (output_dir / "moleculetype.json").write_text(json.dumps(moltype.to_dict()))

    # Read from the file
    new_moleculetype = MoleculeType.from_dict(json.loads((output_dir / "moleculetype.json").read_text()))

    assert moltype.name == new_moleculetype.name
    assert moltype.nrexcl == new_moleculetype.nrexcl
    
    