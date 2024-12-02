#!/usr/bin/env python 
import numpy as np
import pandas as pd
from tqdm import tqdm

import MDAnalysis as mda
from MDAnalysis.analysis import distances
from math import degrees

class Monomer:
    """
    Docstring go brrr
    """
    def __init__(self, monomerName) -> None:
        """
        Create a Monomer MDAnalysis Universe from the filepath to a PDB, and
        save a copy of the Monomer's residues and atoms for quicker access.

        Args:
            monomerName (str): polymer constructed with the Polymer
                    class, supplied from a Polymer.polymer attribute

        Remarks:
            Polymer.dihedral_solver(), Polymer.dist(), Polymer.shuffle() and Polymer.shuffler() rely on connectivity information.
            If your input file does not contain connectivity information (e.g. CONECT records in a pdb), these may not work as intended 
        """
        self.monomer = mda.Universe(monomerName)
        self.residues = self.monomer.residues
        self.atoms = self.monomer.atoms
    
    def select_atoms(self, selection):
        """
        Selection method that selects atoms from the monomer's MDAnalysis
        Universe by leveraging the MDAnalysis atom selection.

        Args:
            selection (str): MDAnalysis atom selection string, for more 
                    information on supported atom selection formats see 
                    [MDAnalysis Atom selection language](https://userguide.mdanalysis.org/stable/selections.html)
        """
        return self.monomer.select_atoms(selection)

