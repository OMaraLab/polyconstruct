import json
from pathlib import Path
from polytop.Angles import Angle
from polytop.Atoms import Atom
from polytop.Dihedrals import Dihedral, Dihedral_type
from polytop.Bonds import Bond
from polytop.Topology import Topology


def test_proper_dihedral_creation()->None:
    atom_a = Atom(
        atom_id=2,
        atom_type="C",
        residue_id=1,
        residue_name="ARG",
        atom_name="C1",
        charge_group_num=2,
        partial_charge=0.5,
        mass=12.011,
    )
    atom_b = Atom(
        atom_id=4,
        atom_type="C",
        residue_id=1,
        residue_name="ARG",
        atom_name="C2",
        charge_group_num=4,
        partial_charge=0.5,
        mass=12.011,
    )
    atom_c = Atom(
        atom_id=8,
        atom_type="C",
        residue_id=1,
        residue_name="ARG",
        atom_name="C3",
        charge_group_num=8,
        partial_charge=0.5,
        mass=12.011,
    )
    atom_d = Atom(
        atom_id=9,
        atom_type="C",
        residue_id=1,
        residue_name="ARG",
        atom_name="C4",
        charge_group_num=9,
        partial_charge=0.5,
        mass=12.011,
    )

    Bond(atom_a, atom_b, bond_type="1", bond_length=1.5, force_constant=1000)
    Bond(atom_b, atom_c, bond_type="1", bond_length=1.5, force_constant=1000)
    Bond(atom_c, atom_d, bond_type="1", bond_length=1.5, force_constant=1000)

    Angle(atom_a, atom_b, atom_c, angle_type="1", angle_value=120.0, force_constant=200)
    Angle(atom_b, atom_c, atom_d, angle_type="1", angle_value=120.0, force_constant=200)

    dihedral = Dihedral(
        atom_a,
        atom_b,
        atom_c,
        atom_d,
        dihedral_type=1,
        phase_angle=180.0,
        force_constant=25.10400,
        multiplicity=3,
    )

    assert dihedral.atom_a.atom_id == 2
    assert dihedral.atom_b.atom_id == 4
    assert dihedral.atom_c.atom_id == 8
    assert dihedral.atom_d.atom_id == 9
    assert dihedral.dihedral_type == Dihedral_type.proper
    assert dihedral.phase_angle == 180.0
    assert dihedral.force_constant == 25.10400
    assert dihedral.multiplicity == 3

    # test contains other atom 
    assert dihedral.contains_atom(atom_a) == True
    assert dihedral.contains_atom(atom_b) == True
    assert dihedral.contains_atom(atom_c) == True
    assert dihedral.contains_atom(atom_d) == True

    # test from atoms
    assert Dihedral.from_atoms(atom_a,atom_b, atom_c, atom_d) == dihedral
    atoms=[atom_a,atom_b,atom_c, atom_d]
    #test serialization works and doesn't create new Dihedrals when there is an existing dihedral
    assert Dihedral.from_dict(dihedral.to_dict(),atoms) == dihedral 

def test_improper_dihedrals_creation()->None:
    atom_a = Atom(
        atom_id=1,
        atom_type="C",
        residue_id=1,
        residue_name="ARG",
        atom_name="C1",
        charge_group_num=1,
        partial_charge=0.5,
        mass=12.011,
    )
    atom_b = Atom(
        atom_id=2,
        atom_type="C",
        residue_id=1,
        residue_name="ARG",
        atom_name="C2",
        charge_group_num=2,
        partial_charge=0.5,
        mass=12.011,
    )
    atom_c = Atom(
        atom_id=3,
        atom_type="C",
        residue_id=1,
        residue_name="ARG",
        atom_name="C3",
        charge_group_num=3,
        partial_charge=0.5,
        mass=12.011,
    )
    atom_d = Atom(
        atom_id=4,
        atom_type="C",
        residue_id=1,
        residue_name="ARG",
        atom_name="C4",
        charge_group_num=4,
        partial_charge=0.5,
        mass=12.011,
    )
    
    Bond(atom_a, atom_b, bond_type="1", bond_length=1.5, force_constant=1000)
    Bond(atom_a, atom_c, bond_type="1", bond_length=1.5, force_constant=1000)
    Bond(atom_a, atom_d, bond_type="1", bond_length=1.5, force_constant=1000)

    Angle(atom_b, atom_a, atom_c, angle_type="1", angle_value=120.0, force_constant=200)
    Angle(atom_b, atom_a, atom_d, angle_type="1", angle_value=120.0, force_constant=200)

    dihedral = Dihedral(
        atom_a,
        atom_b,
        atom_c,
        atom_d,
        dihedral_type=Dihedral_type.improper,
        phase_angle=180.0,
        force_constant=10.5,
        multiplicity=2,
    )

    assert dihedral.atom_a.atom_id == 1
    assert dihedral.atom_b.atom_id == 2
    assert dihedral.atom_c.atom_id == 3
    assert dihedral.atom_d.atom_id == 4
    assert dihedral.dihedral_type == Dihedral_type.improper
    assert dihedral.phase_angle == 180.0
    assert dihedral.force_constant == 10.5
    assert dihedral.multiplicity == 2


def test_dihedrals_storage()->None:
    atom_a = Atom(
        atom_id=1,
        atom_type="C",
        residue_id=1,
        residue_name="ARG",
        atom_name="C1",
        charge_group_num=1,
        partial_charge=0.5,
        mass=12.011,
    )
    atom_b = Atom(
        atom_id=2,
        atom_type="C",
        residue_id=1,
        residue_name="ARG",
        atom_name="C2",
        charge_group_num=2,
        partial_charge=0.5,
        mass=12.011,
    )
    atom_c = Atom(
        atom_id=3,
        atom_type="C",
        residue_id=1,
        residue_name="ARG",
        atom_name="C3",
        charge_group_num=3,
        partial_charge=0.5,
        mass=12.011,
    )
    atom_d = Atom(
        atom_id=4,
        atom_type="C",
        residue_id=1,
        residue_name="ARG",
        atom_name="C4",
        charge_group_num=4,
        partial_charge=0.5,
        mass=12.011,
    )

    topology = Topology()
    topology.add_atom(atom_a)
    topology.add_atom(atom_b)
    topology.add_atom(atom_c)
    topology.add_atom(atom_d)

    bond_ab = Bond(atom_a, atom_b, bond_type="1", bond_length=1.5, force_constant=1000)
    bond_bc = Bond(atom_b, atom_c, bond_type="1", bond_length=1.5, force_constant=1000)
    bond_cd = Bond(atom_c, atom_d, bond_type="1", bond_length=1.5, force_constant=1000)
    
    angle_abc = Angle(
        atom_a, atom_b, atom_c, angle_type="1", angle_value=120.0, force_constant=200
    )
    angle_bcd = Angle(
        atom_b, atom_c, atom_d, angle_type="1", angle_value=120.0, force_constant=200
    )

    dihedral_abcd = Dihedral(
        atom_a,
        atom_b,
        atom_c,
        atom_d,
        dihedral_type=Dihedral_type.proper,
        phase_angle=180.0,
        force_constant=10.5,
        multiplicity=2,
    )

    assert topology.get_bond(1, 2) is not None
    assert topology.get_bond(2, 3) is not None
    assert topology.get_bond(3, 4) is not None
    assert topology.get_angle(1, 2, 3) is not None
    assert topology.get_angle(2, 3, 4) is not None
    dihedral = topology.get_dihedral(1, 2, 3, 4)
    assert dihedral is not None
    assert dihedral.dihedral_type == Dihedral_type.proper

    assert angle_abc in bond_ab.angles
    assert angle_abc in bond_bc.angles
    assert angle_bcd in bond_bc.angles
    assert angle_bcd in bond_cd.angles

    assert dihedral_abcd in angle_abc.dihedrals
    assert dihedral_abcd in angle_bcd.dihedrals

def test_dihedral_serialization(output_dir: Path):
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
            atom_id=2,
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
        ),
        Atom(
            atom_id=4,
            atom_type="C",
            residue_id=1,
            residue_name="ARG",
            atom_name="C9",
            charge_group_num=9,
            partial_charge=0.5,
            mass=12.011,
        )
    ]

    Bond(atomlist[0], atomlist[1], bond_type="1", bond_length=0.147, force_constant=265265.0)
    Bond(atomlist[1], atomlist[2], bond_type="1", bond_length=0.147, force_constant=265265.0)
    Bond(atomlist[2], atomlist[3], bond_type="1", bond_length=0.147, force_constant=265265.0)

    Angle(atomlist[0], atomlist[1], atomlist[2], angle_type="1", angle_value=109.5, force_constant=29288.0)
    Angle(atomlist[1], atomlist[2], atomlist[3], angle_type="1", angle_value=109.5, force_constant=29288.0)
    
    dihedral = Dihedral(atomlist[0], atomlist[1], atomlist[2], atomlist[3], dihedral_type=Dihedral_type.proper, phase_angle=0.0, force_constant=0.0, multiplicity=1)
        
    # Write to the file
    (output_dir / "dihedral.json").write_text(json.dumps(dihedral.to_dict()))

    # Read from the file
    new_dihedral = Dihedral.from_dict(json.loads((output_dir / "dihedral.json").read_text()),atoms=atomlist)
            
    assert dihedral.atom_a.atom_id == new_dihedral.atom_a.atom_id
    assert dihedral.atom_b.atom_id == new_dihedral.atom_b.atom_id
    assert dihedral.atom_c.atom_id == new_dihedral.atom_c.atom_id
    assert dihedral.atom_d.atom_id == new_dihedral.atom_d.atom_id
    assert dihedral.dihedral_type == new_dihedral.dihedral_type
    assert dihedral.phase_angle == new_dihedral.phase_angle
    assert dihedral.force_constant == new_dihedral.force_constant
    assert dihedral.multiplicity == new_dihedral.multiplicity
    
