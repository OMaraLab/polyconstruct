import json
from pathlib import Path
from polytop import *
from polytop.topology import Topology
from polytop.visualize import Visualize
import pytest

def test_invalid_file(data_dir: Path):
    with pytest.raises(FileNotFoundError):
        Topology.from_ITP(data_dir/"non_existent.itp")


def test_invalid_section_warning(data_dir: Path):
    with pytest.warns(UserWarning, match="Unknown section"):
        Topology.from_ITP(data_dir/"invalid_section.itp")


def test_arginine_bonds_angles_dihedrals(data_dir: Path):
    arginine = Topology.from_ITP(data_dir/"arginine.itp")
    assert len(arginine.bonds) == 25
    assert len(arginine.angles) == 41
    assert len(arginine.dihedrals) == 12


def test_glutamine_bonds_angles_dihedrals(data_dir: Path):
    glutamine = Topology.from_ITP(data_dir/"glutamine.itp")
    assert len(glutamine.bonds) == 19
    assert len(glutamine.angles) == 31
    assert len(glutamine.dihedrals) == 10


def test_arginine_molecule_type(data_dir: Path):
    arginine = Topology.from_ITP(data_dir/"arginine.itp")
    assert arginine.molecule_type.name == "AE97"
    assert arginine.molecule_type.nrexcl == 3


def test_glutamine_molecule_type(data_dir: Path):
    glutamine = Topology.from_ITP(data_dir/"glutamine.itp")
    assert glutamine.molecule_type.name == "9NF5"
    assert glutamine.molecule_type.nrexcl == 3


def test_arginine(data_dir: Path, output_dir: Path):
    arginine = Topology.from_ITP(data_dir/"arginine.itp")
    arginine.to_ITP(output_dir/"arginine_new.itp")
    loaded_arginine = Topology.from_ITP(output_dir/"arginine_new.itp")
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


def test_glutamine(data_dir: Path, output_dir: Path):
    glutamine = Topology.from_ITP(data_dir/"glutamine.itp")
    glutamine.to_ITP(output_dir/"glutamine_new.itp")
    loaded_glutamine = Topology.from_ITP(output_dir/"glutamine_new.itp")
    assert len(loaded_glutamine.atoms) == 20
    assert loaded_glutamine.get_atom(20)
    assert loaded_glutamine.get_bond(5, 10)
    assert loaded_glutamine.get_angle(5, 10, 13)

def test_topology_max_atom_index(data_dir: Path, output_dir: Path):
    arginine = Topology.from_ITP(data_dir/"arginine.itp")
    arg_max_atom_index = arginine.max_atom_index()
    assert arg_max_atom_index["H"] == 26
    assert arg_max_atom_index["C"] == 12
    assert arg_max_atom_index["N"] == 6
    assert arg_max_atom_index["O"] == 2

def test_topology_reorder_atom_indexes(data_dir: Path, output_dir: Path):
    arginine = Topology.from_ITP(data_dir/"arginine.itp")
    for element, number in { "H": 26, "C": 12, "N": 6, "O": 2 }.items(): 
        arginine.reorder_atom_indexes(element, number)
    arg_max_atom_index = arginine.max_atom_index()
    assert arg_max_atom_index["H"] == 39
    assert arg_max_atom_index["C"] == 17
    assert arg_max_atom_index["N"] == 9
    assert arg_max_atom_index["O"] == 3
    
def test_topology_auto_rename_atoms(data_dir: Path, output_dir: Path):
    arginine = Topology.from_ITP(data_dir/"arginine.itp")
    arginine.auto_rename_atoms()
    arg_max_atom_index = arginine.max_atom_index()
    assert arg_max_atom_index["H"] == 14
    assert arg_max_atom_index["C"] == 6
    assert arg_max_atom_index["N"] == 4
    assert arg_max_atom_index["O"] == 2

def test_reverse_topology(data_dir: Path, output_dir: Path):
    glu = Topology.from_ITP(data_dir/"glutamine.itp")
    reverse_glu = glu.reverse()
    reverse_glu.to_ITP(output_dir/"reverse_glutamine.itp")
    
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
        
def test_vizualization(data_dir: Path, output_dir: Path):
    import py3Dmol 
    from rdkit import Chem
    from rdkit.Chem import AllChem

    arg = Topology.from_ITP(data_dir/"arginine.itp")
    viewer = py3Dmol.view(width=400, height=400)
    Visualize(arg).draw3D(viewer)

def test_infer_bond_order(data_dir: Path, output_dir: Path):
    arg = Topology.from_ITP(data_dir/"arginine.itp")
    arg = Visualize(arg).topology # visualize constructor infers bond order in it's internal instance of the topology 
    for bond in arg.bonds:
        assert bond.order is not None
        if bond is arg.get_bond(23,24):
            assert bond.order == 2
        elif bond is arg.get_bond(4,8):
            assert bond.order == 2
        else:
            assert bond.order == 1
            
def test_remove(data_dir: Path, output_dir: Path):
    arg = Topology.from_ITP(data_dir/"arginine.itp")
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


def test_serialization(data_dir: Path, output_dir: Path):
    arg = Topology.from_ITP(data_dir/"arginine.itp")
    
    # Write to the file
    (output_dir / "arg.json").write_text(json.dumps(arg.to_dict()))

    # Read from the file
    new_dihedral = Topology.from_dict(json.loads((output_dir/"arg.json").read_text()))
    
    assert len(arg.atoms) == len(new_dihedral.atoms)
    assert len(arg.bonds) == len(new_dihedral.bonds)
    assert len(arg.angles) == len(new_dihedral.angles)
    assert len(arg.dihedrals) == len(new_dihedral.dihedrals)
    assert len(arg.preamble) == len(new_dihedral.preamble)
    assert len(arg.pairs) == len(new_dihedral.pairs)
    assert len(arg.exclusions) == len(new_dihedral.exclusions)
    assert arg.molecule_type == new_dihedral.molecule_type

def is_close(actual,expected) -> bool:
    return abs(actual-expected) < 1e-6
    
    
def test_net_charge(data_dir: Path, output_dir: Path):
    arg = Topology.from_ITP(data_dir/"arginine.itp")
    assert is_close(arg.netcharge, 0)
    # change the net charge 
    arg.netcharge = 1
    assert is_close(arg.netcharge, 1)
    # change the net charge to a fraction of a negative value
    arg.netcharge = -0.5
    assert is_close(arg.netcharge, -0.5)
    arg.to_ITP(output_dir/"arginine_negative_charge.itp")
    
    new_arg = Topology.from_ITP(data_dir/"arginine.itp")
    new_arg.proportional_charge_change(-0.5)
    assert is_close(new_arg.netcharge, -0.5)
    assert not is_close(new_arg.atoms[0].partial_charge, arg.atoms[0].partial_charge)
    arg.to_ITP(output_dir/"arginine_negative_charge_proportional.itp")
    
    
    
