#!/usr/bin/env python 
import random
import pytest
from polyconf.polyconf.Monomer import Monomer
import MDAnalysis as mda
from .test_data import test_pdbs # This is a list of test pdb files that we tested for existence in the test_data.py file


@pytest.mark.parametrize("pdb_filename", test_pdbs) # this parameterization will run this test for every element in test_pbs
def test_monomer_initialization(data_dir, pdb_filename: str):
    """
    Test that we can initialize a Monomer object from a pdb file
    """
    monomer = Monomer(data_dir / pdb_filename)
    assert isinstance(monomer.monomer, mda.Universe), "Monomer.monomer should be an MDAnalysis Universe object"
    # count the lines in the pdb file that begin with HETATM
    with open(data_dir / pdb_filename) as f:
        hetatm_lines = [line for line in f if line.startswith("HETATM")]
    assert len(monomer.atoms) == len(hetatm_lines), "The number of atoms in the Monomer should match the number of HETATM lines in the pdb file"

@pytest.mark.parametrize("pdb_filename", test_pdbs)
def test_select_atoms(data_dir, pdb_filename: str):
    """
    Test that we can find an atom by name in the monomer using the select_atoms method
    """
    monomer = Monomer(data_dir / pdb_filename)
    random_atom = random.choice(monomer.atoms)
    selected_atoms = monomer.select_atoms(f"name {random_atom.name}")
    assert len(selected_atoms) == 1, f"There should only be one {random_atom.name} atom in {pdb_filename}"
    assert selected_atoms[0].name == random_atom.name, f"The selected atom should be named {random_atom.name}"

@pytest.mark.parametrize("pdb_filename", test_pdbs)
def test_monomer_from_u(data_dir, pdb_filename: str):
    """
    Test that the monomer_from_u Monomer Class method works
    """
    universe = mda.Universe(data_dir / pdb_filename)
    monomer = Monomer.monomer_from_u(universe)
    assert len(monomer.atoms) == len(universe.atoms)