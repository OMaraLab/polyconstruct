import pytest

test_pdbs = ['CODI_bonds.pdb']

@pytest.mark.parametrize("pdb_filename", test_pdbs)
def test_data_files_exist(data_dir, pdb_filename):
    assert (data_dir / pdb_filename).exists(), f"{pdb_filename} does not exist in {data_dir}"
