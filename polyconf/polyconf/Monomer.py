#!/usr/bin/env python 
import MDAnalysis as mda

class Monomer:
    """
    Create a Monomer MDAnalysis Universe from the provided PDB filepath.

    Remarks:
            :func:`Polymer.dihedral_solver`, :func:`Polymer.dist`, 
            :func:`Polymer.shuffle` and :func:`Polymer.shuffler` rely on 
            connectivity information. If your input file does not contain
            connectivity information (e.g. CONECT records in a pdb),
            these may not work as intended.
    """
    def __init__(self, monomerName:str) -> None:
        """
        Create a Monomer MDAnalysis Universe from the filepath to a PDB, and
        saves a copy of the Monomer's residues and atoms for quicker access.
        
        :param monomerName: filepath to the monomer .pdb file
        :type monomerName: str
        """
        self.monomer = mda.Universe(monomerName)
        self.residues = self.monomer.residues
        self.atoms = self.monomer.atoms
    
    def select_atoms(self, selection:str) -> mda.AtomGroup:
        """
        Selection method that selects atoms from the monomer's MDAnalysis
        Universe by leveraging the MDAnalysis atom selection.
        MDAnalysis atom selection string, for more
        information on supported atom selection formats see

        :param selection: MDAnalysis atom selection string, for more 
            information on supported atom selection formats see 
            `MDAnalysis Atom selection language <https://userguide.mdanalysis.org/stable/selections.html>`
        :type selection: str
        :return: An MDAnalysis AtomGroup containing only the selected atoms
        :rtype: mda.AtomGroup
        """
        return self.monomer.select_atoms(selection)