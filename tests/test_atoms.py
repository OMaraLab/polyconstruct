import json
from polytop.atoms import Atom
from polytop.topology import Topology


def test_atom_creation()->None:
    atom = Atom(
        atom_id=26,
        atom_type="N3",
        residue_id=1,
        residue_name="ARG",
        atom_name="NH2",
        charge_group_num=26,
        partial_charge=-0.940000,
        mass=14.0067,
    )
    assert atom.atom_id == 26
    assert atom.atom_type == "N3"
    assert atom.residue_id == 1
    assert atom.residue_name == "ARG"
    assert atom.atom_name == "NH2"
    assert atom.charge_group_num == 26
    assert atom.partial_charge == -0.940000
    assert atom.mass == 14.0067

def test_repr():
    arg=Topology.from_ITP("tests/samples/arginine.itp")
    assert arg.atoms[0].__repr__()=='Atom(1H)->[2] charge=0.412'
    
def test_pseudoatoms():
    atom = Atom(
        atom_id=1,
        atom_type="X",
        residue_id=1,
        residue_name="ARG",
        atom_name="X1",
        charge_group_num=26,
        partial_charge=-0.940000,
        mass=14.0067,
    )
    assert atom.is_virtual
    
def test_atom_creation_from_line():
    atom = Atom.from_line("    1  HS14    1    AE97    H26    1    0.412   1.0080")
    assert atom.atom_id == 1
    assert atom.atom_type == "HS14"
    assert atom.residue_id == 1
    assert atom.residue_name == "AE97"
    assert atom.atom_name == "H26"
    assert atom.charge_group_num == 1
    assert atom.partial_charge == 0.412
    assert atom.mass == 1.0080
    
def test_atom_serialization():
    atom = Atom(
        atom_id=26,
        atom_type="N3",
        residue_id=1,
        residue_name="ARG",
        atom_name="NH2",
        charge_group_num=26,
        partial_charge=-0.940000,
        mass=14.0067,
    )
    atom_dict = atom.to_dict()
    with open("tests/samples/atom.json", "w") as f:
        json.dump(atom_dict, f)
    with open("tests/samples/atom.json", "r") as f:
        new_atom = Atom.from_dict(json.load(f))
    assert atom.atom_id == new_atom.atom_id
    assert atom.atom_type == new_atom.atom_type
    assert atom.residue_id == new_atom.residue_id
    assert atom.residue_name == new_atom.residue_name
    assert atom.atom_name == new_atom.atom_name
    assert atom.charge_group_num == new_atom.charge_group_num
    assert atom.partial_charge == new_atom.partial_charge
    assert atom.mass == new_atom.mass
    