import random
import pytest
import MDAnalysis as mda
from polyconf.PDB import PDB
from .test_data import test_pdbs

@pytest.mark.parametrize("pdb_file", test_pdbs)
def test_pdb_initialization(data_dir, pdb_file):
    polymer = mda.Universe(f"{data_dir}/{pdb_file}")
    pdb = PDB(polymer)
    assert pdb.polymer == polymer
    assert len(pdb.atoms) == len(polymer.atoms)

@pytest.mark.parametrize("pdb_file", test_pdbs)
def test_cleanup(data_dir, pdb_file):
    polymer = mda.Universe(f"{data_dir}/{pdb_file}")
    pdb = PDB(polymer)
    pdb.cleanup()
    box = (pdb.atoms.positions).max(axis=0) - (pdb.atoms.positions).min(axis=0) + [10, 10, 10]
    assert pdb.dimensions[:3] == pytest.approx(box, rel=1e-1), "Box is not at least 10 Angstroms larger than the polymer"
    assert pdb.dimensions[3:] == [90, 90, 90], "Box angles are not 90 degrees"
    assert (pdb.atoms.positions >= -box/2).all(), "Some atoms are located outside the box"

@pytest.mark.parametrize("pdb_file", test_pdbs)
def test_write(data_dir, pdb_file, tmp_path):
    polymer = mda.Universe(f"{data_dir}/{pdb_file}")
    pdb = PDB(polymer)
    output_file = tmp_path / "output.pdb"
    for atom in pdb.polymer.atoms: # Test for every atom in the polymer
        pdb._write(f"name {atom.name}", str(output_file))
        assert output_file.exists()
        with open(output_file, 'r') as f:
            content = f.read()
        assert f" {atom.name} " in content, f"Atom {atom.name} not found in output file"

@pytest.mark.parametrize("pdb_file", test_pdbs)
def test_save(data_dir, pdb_file, tmp_path):
    polymer = mda.Universe(f"{data_dir}/{pdb_file}")
    pdb = PDB(polymer)
    output_file = tmp_path / "output.gro"
    random_atom = random.choice(pdb.polymer.atoms) # just test one atom at random 
    pdb.save(dummyAtoms=random_atom.name, fname=str(output_file), selectionString="all", pdb=False)
    assert output_file.exists()
    with open(output_file, 'r') as f:
        content = f.read()
    assert f" {random_atom.name} " in content, f"Atom {random_atom.name} not found in output file"
