#!/usr/bin/env python 
import numpy as np
import pandas as pd
from tqdm import tqdm  # progress bar

import MDAnalysis as mda
from MDAnalysis.analysis import distances
from math import degrees

import random # we'll be using this to randomise monomer positions

import networkx as nx  # we'll be using this to define halves of the molecule during dihedral shuffling

class Monomer:
    """
    Docstring go brrr
    """
    def __init__(self, monomerName) -> None:
        self.monomer = mda.Universe(monomerName)
        self.residues = self.monomer.residues
        self.atoms = self.monomer.atoms
    
    def select_atoms(self, selection):
        return self.monomer.select_atoms(selection)

