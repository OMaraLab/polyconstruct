#!/usr/bin/python
import sys
MIN_PYTHON = (3, 8)
if sys.version_info < MIN_PYTHON:
    sys.exit("Python %s.%s or later is required.\n" % MIN_PYTHON)

from .molecule import Molecule
from .monomer import Monomer
from .polymer import Polymer

def polytopVer():
    return 0.1

