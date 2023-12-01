def test_data_files_exist(data_dir):
    assert (data_dir / 'arginine.itp').exists(), "arginine.itp does not exist"
    assert (data_dir / 'glutamine.itp').exists(), "glutamine.itp does not exist"
