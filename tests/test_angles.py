import json
from polytop.angles import Angle
from polytop.atoms import Atom
from polytop.bonds import Bond


def test_angle_creation()->None:
    atom_a = Atom(
        atom_id=9,
        atom_type="C",
        residue_id=1,
        residue_name="ARG",
        atom_name="C9",
        charge_group_num=9,
        partial_charge=0.5,
        mass=12.011,
    )
    atom_b = Atom(
        atom_id=12,
        atom_type="C",
        residue_id=1,
        residue_name="ARG",
        atom_name="C12",
        charge_group_num=12,
        partial_charge=0.5,
        mass=12.011,
    )
    atom_c = Atom(
        atom_id=15,
        atom_type="C",
        residue_id=1,
        residue_name="ARG",
        atom_name="C15",
        charge_group_num=15,
        partial_charge=0.5,
        mass=12.011,
    )
    bond_ab = Bond(atom_a, atom_b, bond_type="1", bond_length=1.5, force_constant=1000)
    bond_bc = Bond(atom_b, atom_c, bond_type="1", bond_length=1.5, force_constant=1000)
    
    angle = Angle(
        atom_a,
        atom_b,
        atom_c,
        angle_type="1",
        angle_value=109.5,
        force_constant=29288.0,
    )
    assert angle.atom_a.atom_id == 9
    assert angle.atom_b.atom_id == 12
    assert angle.atom_c.atom_id == 15
    assert angle.angle_type == "1"
    assert angle.angle_value == 109.5
    assert angle.force_constant == 29288.0

def test_angle_serialization():
    atomlist = [
        Atom(
            atom_id=1,
            atom_type="C",
            residue_id=1,
            residue_name="ARG",
            atom_name="C9",
            charge_group_num=9,
            partial_charge=0.5,
            mass=12.011,
        ),
        Atom(
            atom_id=4,
            atom_type="C",
            residue_id=1,
            residue_name="ARG",
            atom_name="C12",
            charge_group_num=12,
            partial_charge=0.5,
            mass=12.011,
        ),
        Atom(
            atom_id=3,
            atom_type="C",
            residue_id=1,
            residue_name="ARG",
            atom_name="C15",
            charge_group_num=15,
            partial_charge=0.5,
            mass=12.011,
        )
    ]

    Bond(atomlist[0], atomlist[1], bond_type="1", bond_length=0.147, force_constant=265265.0)
    Bond(atomlist[1], atomlist[2], bond_type="1", bond_length=0.147, force_constant=265265.0)

    angle = Angle(atomlist[0], atomlist[1], atomlist[2], angle_type="1", angle_value=109.5, force_constant=29288.0)
    angle_dict = angle.to_dict()
    with open("tests/samples/angle.json", "w") as f:
        json.dump(angle_dict, f)
    with open("tests/samples/angle.json", "r") as f:
        new_bond = Angle.from_dict(json.load(f),atoms=atomlist)
    assert angle.atom_a.atom_id == new_bond.atom_a.atom_id
    assert angle.atom_b.atom_id == new_bond.atom_b.atom_id
    assert angle.atom_c.atom_id == new_bond.atom_c.atom_id
    assert angle.angle_type == new_bond.angle_type
    assert angle.angle_value == new_bond.angle_value
    assert angle.force_constant == new_bond.force_constant
    