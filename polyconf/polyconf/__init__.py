#!/usr/bin/python
"""
PolyConf: A Python library for creating polymers from monomers using GROMACS PDB files

Copyright (C) 2024 O'Mara Group
License: GPL
Website: https://aibn.uq.edu.au/omara
Author: Ada Quinn and Luna Morrow
"""

__version__ = "1.0"
__author__ = "Ada Quinn and Luna Morrow"
__license__ = "GPL"
__website__ = "https://aibn.uq.edu.au/omara"

from .Monomer import Monomer
from .Polymer import Polymer
from .PDB import PDB
from .polyconf_automatic import main