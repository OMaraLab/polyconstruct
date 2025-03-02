#!/usr/bin/env python 
from __future__ import annotations
import MDAnalysis as mda
from MDAnalysis import Merge

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
    def __init__(self, monomerName, fromUni = False) -> None:
        """
        Create a Monomer MDAnalysis Universe from the filepath to a PDB, and
        saves a copy of the Monomer's residues and atoms for quicker access.
        
        :param monomerName: filepath to the monomer .pdb file
        :type monomerName: str
        """
        if fromUni:
            self.monomer = monomerName
        else:
            self.monomer = mda.Universe(monomerName)
        self.residues = self.monomer.residues
        self.atoms = self.monomer.atoms
        self.dimensions = self.monomer.dimensions

    @classmethod
    def monomer_from_u(cls, universe:mda.Universe) -> 'Monomer':
        """
        Create Monomer from an MDAnalysis Universe

        :param universe: _description_
        :type universe: mda.Universe
        :return: _description_
        :rtype: Monomer
        """
        return cls(universe.select_atoms("all"), fromUni=True)
    
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

    def copy(self) -> Monomer:
        """
        Use MDAnalysis Universe.copy() to create a new MDAnalysis Universe that
        is an exact copy of this one containing the polymer. Deepcopy is not
        used as it can be problematic with open file sockets.

        :return: A new Polymer that is an exact copy of this polymer
        :rtype: Polymer
        """
        new = self.monomer.copy()
        return Monomer.monomer_from_u(new)

    # def __eq__(self, value):
    #     new = self.atoms.subtract(value.atoms)
    #     print(len(new))
    #     return self.select_atoms("all") == value.select_atoms("all") # and self.residues.equals(value.residues)