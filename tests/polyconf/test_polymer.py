#!/usr/bin/env python 
import pytest
import random
from polyconf.polyconf.Polymer import Polymer
from polyconf.polyconf.Monomer import Monomer
import MDAnalysis as mda
from .test_data import test_pdbs

@pytest.mark.parametrize("pdb_file", test_pdbs)
def test_polymer_initialization(data_dir, pdb_file):
    '''
    Test that we can initialize a Polymer object from a pdb file
    And that it defines a bunch of attributes and functions that we expect
    '''
    
    first_monomer = Monomer(f"{data_dir}/{pdb_file}")
    polymer = Polymer(first_monomer)
    assert polymer.first == first_monomer
    assert polymer.first.residues.resids[0] == 1

    # let's make sure some attributes we expect are set
    assert hasattr(polymer, 'first')
    assert hasattr(polymer, 'polymer')
    assert hasattr(polymer, 'atoms')

    # let's make sure some functions we expect are defined
    methods = ['select_atoms', 'copy', 'renamer', 'extend', 'dist', 'gencomp',
               'shuffle', 'shuffler', 'dihedral_solver', 'gen_pairlist',
               'newresid', 'maxresid']
    for method in methods:
        assert hasattr(Polymer, method), f"Polymer class should have a method named {method}"

def test_polymer_initialization_arguments():
    '''
    Test that the Polymer class raises an error when initialized with a non-Monomer object
    '''
    with pytest.raises(ValueError):
        Polymer(None)
    with pytest.raises(TypeError):
        Polymer("not a Monomer")

@pytest.mark.parametrize("pdb_file", test_pdbs)
def test_polymer_attributes(data_dir, pdb_file):
    first_monomer = Monomer(f"{data_dir}/{pdb_file}")
    polymer = Polymer(first_monomer)
    assert polymer.first.residues.resids[0] == 1

@pytest.mark.parametrize("pdb_file", test_pdbs)
def test_select_atoms(data_dir, pdb_file: str):
    """
    Test that we can find an atom by name in the polymer using the select_atoms method
    """
    polymer = Polymer(Monomer(data_dir / pdb_file))
    random_atom = random.choice(polymer.atoms)
    selected_atoms = polymer.select_atoms(f"name {random_atom.name}")
    assert len(selected_atoms) == 1, f"There should only be one {random_atom.name} atom in {pdb_file}"
    assert selected_atoms[0].name == random_atom.name, f"The selected atom should be named {random_atom.name}"

# Add more tests for other methods and functionalities of the Polymer class as needed