from pathlib import Path
from polytop import *
import subprocess
import os

import pytest

@pytest.mark.xfail(reason="Requires a docker container with Gromacs installed")
def test_gromacs_help(data_dir, output_dir):
    simulation = Gromacs(data_dir, output_dir)
    simulation.run(["gmx", "help", "commands"], "help.txt")

    # Assert that the output file was created
    output_file = output_dir / 'help.txt'
    assert output_file.exists(), f"{output_file} does not exist"
        