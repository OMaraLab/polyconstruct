# hand polymer ITP connectivity info to polyconf, instead of relying on input pdb's having connectivity

#!/usr/bin/env python 
import numpy as np
import pandas as pd
from tqdm import tqdm  # progress bar

import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.core.groups import AtomGroup
from math import degrees

import random # we'll be using this to randomise monomer positions

import networkx as nx  # we'll be using this to define halves of the molecule during dihedral shuffling

class PDB:
    """
    docstring go brr
    """
    def __init__(self, polymer) -> None:
        self.polymer = polymer
        self.atoms = polymer.atoms
    
    def cleanup(self):
        # adjust box size, center polymer in box
        box = (self.atoms.positions).max(axis=0)- (self.atoms.positions).min(axis=0) + [10,10,10]
        self.dimensions = list(box)  + [90]*3
        
        cog = self.atoms.center_of_geometry(wrap=False)
        box_center = box / 2
        self.atoms.translate(box_center - cog)

    def write(self, selection, name):
        atoms = self.polymer.select_atoms(selection) #AtomGroup(atoms)
        atoms.write(name)

    def save(self,fname='polymer.gro', selectionString = None):
        if selectionString:
            self.select_atoms(f'{selectionString} and not name CN CMA CP CQ').atoms.write(fname)  # save, excluding dummy atoms
        else:
            # atoms = self.select_atoms('not name CN CMA CP CQ')
            self.write('not name CN CMA CP CQ', fname)  # save, excluding dummy atoms

    def crudesave(self,fname='polymer_crude.gro'):
        # save including dummy atoms, useful for debugging the geometry transforms
        self.atoms.write(fname)

    def select_atoms(self, selection):
        return self.polymer.select_atoms(selection)
