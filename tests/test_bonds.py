import json
from polytop.atoms import Atom
from polytop.bonds import Bond
from polytop.topology import Topology


def test_bond_creation()->None:
    atom_a = Atom(
        atom_id=20,
        atom_type="C",
        residue_id=1,
        residue_name="ARG",
        atom_name="C20",
        charge_group_num=20,
        partial_charge=0.5,
        mass=12.011,
    )
    atom_b = Atom(
        atom_id=22,
        atom_type="C",
        residue_id=1,
        residue_name="ARG",
        atom_name="C22",
        charge_group_num=22,
        partial_charge=0.5,
        mass=12.011,
    )
    bond = Bond(
        atom_a, atom_b, bond_type="1", bond_length=0.147, force_constant=265265.0
    )
    assert bond.atom_a.atom_id == 20
    assert bond.atom_b.atom_id == 22
    assert bond.bond_type == "1"
    assert bond.bond_length == 0.147
    assert bond.force_constant == 265265.0

def test_bond_traversal():
    arg = Topology.from_ITP("tests/samples/arginine.itp")
    bond = arg.get_bond(12, 15)
    assert bond.atom_a.atom_id == 12
    assert bond.atom_b.atom_id == 15
    LHS = bond.LHS()
    RHS = bond.RHS()
    assert len(LHS) == 14
    assert len(RHS) == 12
    allatoms = set(arg.atoms)
    assert allatoms == LHS | RHS
    assert LHS & RHS == set()

def test_bond_serialization():
    atomlist = [
        Atom(
            atom_id=1,
            atom_type="N3",
            residue_id=1,
            residue_name="test",
            atom_name="NH2",
            charge_group_num=26,
            partial_charge=-0.940000,
            mass=14.0067,
        ),
        Atom(
            atom_id=2,
            atom_type="H1",
            residue_id=1,
            residue_name="test",
            atom_name="NH2",
            charge_group_num=26,
            partial_charge=-0.940000,
            mass=14.0067,
        )
    ]
    bond = Bond(atomlist[0], atomlist[1], bond_type="1", bond_length=0.147, force_constant=265265.0)
    bond_dict = bond.to_dict()
    with open("tests/samples/bonds.json", "w") as f:
        json.dump(bond_dict, f)
    with open("tests/samples/bonds.json", "r") as f:
        new_bond = Bond.from_dict(json.load(f),atoms=atomlist)
    assert new_bond.atom_a.atom_id == bond.atom_a.atom_id
    assert new_bond.atom_b.atom_id == bond.atom_b.atom_id
    assert new_bond.bond_type == bond.bond_type
    assert new_bond.bond_length == bond.bond_length
    assert new_bond.force_constant == bond.force_constant

def test_arg_double_bonds():
    arg = Topology.from_ITP("tests/samples/arginine.itp")
    double_bond_1 = arg.get_bond('C11','O2')
    assert double_bond_1.order == 1
    assert double_bond_1
    double_bond_1.order = 2 # ITP files do not include bond orders
    assert double_bond_1.order == 2

    double_bond_2 = arg.get_bond('C12','N4')
    assert double_bond_2.order == 1
    assert double_bond_2
    double_bond_2.order = 2 # ITP files do not include bond orders
    assert double_bond_2.order == 2
