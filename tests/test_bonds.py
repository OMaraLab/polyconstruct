import json
from pathlib import Path

import pytest
from polytop.polytop.Atoms import Atom
from polytop.polytop.Bonds import Bond
from polytop.polytop.Topology import Topology


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
    # test otheratom across from a bond
    assert bond.other_atom(atom_a) == atom_b
    assert bond.other_atom(atom_b) == atom_a
    # test from_atoms
    assert Bond.from_atoms(atom_a,atom_b) == bond 
    atoms=[atom_a,atom_b]
    #test serialization works and doesn't create new Bonds when there is an existing bond
    assert Bond.from_dict(bond.to_dict(),atoms) == bond 

def test_deduplicate_bonds():
    atom1 = Atom(
        atom_id=20,
        atom_type="C",
        residue_id=1,
        residue_name="ARG",
        atom_name="C20",
        charge_group_num=20,
        partial_charge=0.5,
        mass=12.011,
    )
    atom2 = Atom(
        atom_id=22,
        atom_type="C",
        residue_id=1,
        residue_name="ARG",
        atom_name="C22",
        charge_group_num=22,
        partial_charge=0.5,
        mass=12.011,
    )
    bond1 = Bond(
        atom1, atom2, bond_type="1", bond_length=0.147, force_constant=265265.0
    )
    bond2 = Bond(
        atom1, atom2, bond_type="1", bond_length=0.147, force_constant=265265.0
    )
    assert len(atom1.bonds) == 2
    atom1.deduplicate_bonds()
    assert len(atom1.bonds) == 1    

def test_bond_traversal(data_dir: Path):
    arg = Topology.from_ITP(data_dir / "arginine.itp")
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

def test_bond_serialization(output_dir: Path):
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
        
    # Write to the file
    (output_dir / "bonds.json").write_text(json.dumps(bond.to_dict()))

    # Read from the file
    new_bond = Bond.from_dict(json.loads((output_dir / "bonds.json").read_text()),atoms=atomlist)
            
    assert new_bond.atom_a.atom_id == bond.atom_a.atom_id
    assert new_bond.atom_b.atom_id == bond.atom_b.atom_id
    assert new_bond.bond_type == bond.bond_type
    assert new_bond.bond_length == bond.bond_length
    assert new_bond.force_constant == bond.force_constant

def test_arg_double_bonds(data_dir: Path):
    arg = Topology.from_ITP(data_dir / "arginine.itp")
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
