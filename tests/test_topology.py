import json
from polytop import *


import pytest


def test_invalid_file():
    with pytest.raises(FileNotFoundError):
        Topology.from_ITP("tests/samples/non_existent.itp")


def test_invalid_section_warning():
    with pytest.warns(UserWarning, match="Unknown section"):
        Topology.from_ITP("tests/samples/invalid_section.itp")


def test_arginine_bonds_angles_dihedrals():
    arginine = Topology.from_ITP("tests/samples/arginine.itp")
    assert len(arginine.bonds) == 25
    assert len(arginine.angles) == 41
    assert len(arginine.dihedrals) == 12


def test_glutamine_bonds_angles_dihedrals():
    glutamine = Topology.from_ITP("tests/samples/glutamine.itp")
    assert len(glutamine.bonds) == 19
    assert len(glutamine.angles) == 31
    assert len(glutamine.dihedrals) == 10


def test_arginine_molecule_type():
    arginine = Topology.from_ITP("tests/samples/arginine.itp")
    assert arginine.molecule_type.name == "AE97"
    assert arginine.molecule_type.nrexcl == 3


def test_glutamine_molecule_type():
    glutamine = Topology.from_ITP("tests/samples/glutamine.itp")
    assert glutamine.molecule_type.name == "9NF5"
    assert glutamine.molecule_type.nrexcl == 3


def test_arginine():
    arginine = Topology.from_ITP("tests/samples/arginine.itp")
    arginine.to_ITP("tests/samples/arginine_new.itp")
    loaded_arginine = Topology.from_ITP("tests/samples/arginine_new.itp")
    assert (
        len(loaded_arginine.atoms) == 26
    )  # verify that the number of atoms is correct
    assert loaded_arginine.get_atom(26)  # verify that the atom with id 26 exists
    assert loaded_arginine.get_bond(
        20, 22
    )  # verify that the bond between atoms 20 and 22 exists
    assert loaded_arginine.get_angle(
        9, 12, 15
    )  # verify that the angle between atoms 9, 12 and 15 exists
    assert loaded_arginine.get_dihedral(
        4, 2, 5, 8
    )  # verify that the improper dihedral between atoms 4, 2, 5 and 8 exists
    assert loaded_arginine.get_dihedral(
        2, 4, 8, 9
    )  # verify that the proper dihedral between atoms 2, 4, 8 and 9 exists


def test_glutamine():
    glutamine = Topology.from_ITP("tests/samples/glutamine.itp")
    glutamine.to_ITP("tests/samples/glutamine_new.itp")
    loaded_glutamine = Topology.from_ITP("tests/samples/glutamine_new.itp")
    assert len(loaded_glutamine.atoms) == 20
    assert loaded_glutamine.get_atom(20)
    assert loaded_glutamine.get_bond(5, 10)
    assert loaded_glutamine.get_angle(5, 10, 13)

def test_reverse_topology():
    glu = Topology.from_ITP("tests/samples/glutamine.itp")
    reverse_glu = glu.reverse()
    reverse_glu.to_ITP("tests/samples/reverse_glutamine.itp")
    
    def same_properties(atom_a,atom_b)-> bool: 
        if glu_atom.atom_type != rev_glu_atom.atom_type:
            return False
        if glu_atom.residue_name != rev_glu_atom.residue_name:
            return False
        if glu_atom.residue_id != rev_glu_atom.residue_id:
            return False
        # assert arg_atoms[i].charge_group_num == rev_arg_atom.charge_group_num # it is not dependent on order but grouping
        if glu_atom.partial_charge != rev_glu_atom.partial_charge:
            return False
        if glu_atom.mass != rev_glu_atom.mass:
            return False
        return True
        
    glu_atoms = glu.atoms
    reverse_glu_atoms = reverse_glu.atoms
    assert len(glu_atoms) == len(reverse_glu_atoms)
    for i in range(len(glu_atoms)):
        glu_atom = glu_atoms[i]
        rev_glu_atom = reverse_glu_atoms[-i-1]
        assert same_properties(glu_atom, rev_glu_atom)
        for atom in glu_atom.bond_neighbours():
            assert any(same_properties(atom, rev_atom) for rev_atom in rev_glu_atom.bond_neighbours())
        for atom in glu_atom.angle_neighbours():
            assert any(same_properties(atom, rev_atom) for rev_atom in rev_glu_atom.angle_neighbours())
        for atom in glu_atom.dihedral_neighbours():
            assert any(same_properties(atom, rev_atom) for rev_atom in rev_glu_atom.dihedral_neighbours())
        
def test_vizualization():
    import py3Dmol 
    from rdkit import Chem
    from rdkit.Chem import AllChem

    arg = Topology.from_ITP("tests/samples/arginine.itp")
    viewer = py3Dmol.view(width=400, height=400)
    Visualize(arg).create_py3Dmol_view(viewer)

def test_infer_bond_order():
    arg = Topology.from_ITP("tests/samples/arginine.itp")
    Visualize(arg).infer_bond_orders()
    for bond in arg.bonds:
        assert bond.order is not None
        if bond is arg.get_bond(23,24):
            assert bond.order == 2
        elif bond is arg.get_bond(4,8):
            assert bond.order == 2
        else:
            assert bond.order == 1
            
def test_remove():
    arg = Topology.from_ITP("tests/samples/arginine.itp")
    assert len(arg.atoms) == 26
    assert len(arg.bonds) == 25
    assert len(arg.angles) == 41
    assert len(arg.dihedrals) == 12
    atom = arg.get_atom(1)
    arg.remove_atom(atom)
    assert len(arg.atoms) == 25
    assert len(arg.bonds) == 24
    assert len(arg.angles) == 39
    assert len(arg.dihedrals) == 12
    atom = arg.get_atom(2)
    arg.remove_atom(atom)
    assert len(arg.atoms) == 24
    assert len(arg.bonds) == 22
    assert len(arg.angles) == 36
    assert len(arg.dihedrals) == 9
    atom = arg.get_atom(26)
    arg.remove_atom(atom)
    assert len(arg.atoms) == 23
    assert len(arg.bonds) == 21
    assert len(arg.angles) == 35
    assert len(arg.dihedrals) == 8

def test_split():
    arg = Topology.from_ITP("tests/samples/arginine.itp")
    first_polymerization_bond = arg.get_bond_by_name('N3','H20') # H20 replaced with X)
    LHS,RHS = arg.split(first_polymerization_bond)
    assert len(LHS.atoms) == 26
    assert len(RHS.atoms) == 2
    assert LHS.atoms[21-1].is_virtual
    assert RHS.atoms[0].is_virtual
    arg = Topology.from_ITP("tests/samples/arginine.itp")
    second_polymerization_bond = arg.get_bond_by_name('C11','O1') # O1 replaced with X)
    LHS,RHS = arg.split(second_polymerization_bond)
    assert len(LHS.atoms) == 25
    assert len(RHS.atoms) == 3
    assert LHS.atoms[24].is_virtual
    assert RHS.atoms[0].is_virtual

def test_serialization():
    arg = Topology.from_ITP("tests/samples/arginine.itp")
    with open("tests/samples/arg.json", "w") as f:
        json.dump(arg.to_dict(), f)
    with open("tests/samples/arg.json", "r") as f:
        new_dihedral = Topology.from_dict(json.load(f))
    assert len(arg.atoms) == len(new_dihedral.atoms)
    assert len(arg.bonds) == len(new_dihedral.bonds)
    assert len(arg.angles) == len(new_dihedral.angles)
    assert len(arg.dihedrals) == len(new_dihedral.dihedrals)
    assert len(arg.preamble) == len(new_dihedral.preamble)
    assert len(arg.pairs) == len(new_dihedral.pairs)
    assert len(arg.exclusions) == len(new_dihedral.exclusions)
    assert arg.molecule_type == new_dihedral.molecule_type
