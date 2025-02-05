#!/usr/bin/env python 
import random
import pytest
from polyconf.PDB import PDB
from polyconf.Polymer import Polymer
from polyconf.Monomer import Monomer
import MDAnalysis as mda
from .test_data import test_pdbs

@pytest.mark.parametrize("pdb_file", test_pdbs)
def test_pdb_initialization(data_dir, pdb_file):
    polymer = Polymer(Monomer(f"{data_dir}/{pdb_file}"))
    pdb = PDB(polymer)
    assert pdb.polymer == polymer
    assert len(pdb.atoms) == len(polymer.atoms)

@pytest.mark.parametrize("pdb_file", test_pdbs)
def test_cleanup(data_dir, pdb_file):
    polymer = Polymer(Monomer(f"{data_dir}/{pdb_file}"))
    pdb = PDB(polymer)
    pdb.cleanup()
    box = (pdb.atoms.positions).max(axis=0) - (pdb.atoms.positions).min(axis=0) + [10, 10, 10]
    assert pdb.dimensions[:3] == pytest.approx(box, rel=1e-1), "Box is not at least 10 Angstroms larger than the polymer"
    assert pdb.dimensions[3:] == [90, 90, 90], "Box angles are not 90 degrees"
    assert (pdb.atoms.positions >= -box/2).all(), "Some atoms are located outside the box"

@pytest.mark.parametrize("pdb_file", test_pdbs)
def test_crude_save(data_dir, pdb_file, output_dir):
    polymer = Polymer(Monomer(f"{data_dir}/{pdb_file}"))
    pdb = PDB(polymer)
    output_file = output_dir / "output.pdb"
    pdb.crudesave(fname=str(output_dir / "output"))
    assert output_file.exists()
    with open(output_file, 'r') as f:
        content = f.read()
    for atom in pdb.polymer.atoms: # Test for every atom in the polymer
        assert f" {atom.name} " in content, f"Atom {atom.name} not found in output file"

@pytest.mark.parametrize("pdb_file", test_pdbs)
def test_save(data_dir, pdb_file, output_dir):
    polymer = Polymer(Monomer(f"{data_dir}/{pdb_file}"))
    pdb = PDB(polymer)
    output_file = output_dir / "output.gro"
    pdb.save(fname=str(output_dir / "output"), gmx=True)
    assert output_file.exists()
    with open(output_file, 'r') as f:
        content = f.read()
    for atom in pdb.polymer.atoms: # Test for every atom in the polymer
        assert f" {atom.name} " in content, f"Atom {atom.name} not found in output file"
