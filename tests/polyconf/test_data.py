import pytest

test_pdbs = ['CODI_bonds.pdb']

@pytest.mark.parametrize("pdb_filename", test_pdbs)
def test_data_files_exist(data_dir, pdb_filename):
    assert (data_dir / pdb_filename).exists(), f"{pdb_filename} does not exist in {data_dir}"

PEI_pdbs = ['PEI_start.pdb', 'PEI_end.pdb', 'PEI_monomer.pdb']
PEI_start = ['PEI_start.pdb']
PEI_monomer = ['PEI_monomer.pdb']
PEI_end = ['PEI_end.pdb']

@pytest.mark.parametrize("PEI_filename", PEI_pdbs)
def test_PEI_files_exist(data_dir, PEI_filename):
    assert (data_dir / PEI_filename).exists(), f"{PEI_filename} does not exist in {data_dir}"