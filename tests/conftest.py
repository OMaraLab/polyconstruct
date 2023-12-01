import pytest
from pathlib import Path

@pytest.fixture(scope="session")
def data_dir():
    """Return the path to the data directory"""
    return Path(__file__).parent / 'data'

@pytest.fixture(scope="session")
def output_dir():
    """Common path to the test output directory (not )"""
    output_dir = Path(__file__).parent / 'output'
    output_dir.mkdir(exist_ok=True)
    return output_dir