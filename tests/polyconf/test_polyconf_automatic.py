import pytest
import subprocess
from pathlib import Path


def test_polyconf_automatic_help():
    '''
    Test that the cli script for polyconf_automatic.py runs and prints help message
    '''
    result = subprocess.run(['polyconf-automatic', '--help'], capture_output=True, text=True)
    assert result.returncode == 0
    assert "usage:" in result.stdout
