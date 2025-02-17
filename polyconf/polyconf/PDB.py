#!/usr/bin/env python

import MDAnalysis as mda
    
class PDB:
    """
    Class containing methods to copy and save a :class:`Polymer` to a PDB or
    GMX file.
    """
    def __init__(self, polymer) -> None:
        """
        Save a copy of the Polymer.polymer MDAnalysis Universe and its atoms
        for cleanup and PDB export.

        :param polymer: polymer constructed with :class:`Polymer`, supplied
                from a Polymer.polymer attribute
        :type polymer: :class:`Polymer`
        """
        self.polymer = polymer
        self.atoms = polymer.atoms
    
    def cleanup(self):
        """
        Adjust the box size and center the polymer in the box
        """
        box = (self.atoms.positions).max(axis=0)- (self.atoms.positions).min(axis=0) + [10,10,10]
        self.dimensions = list(box)  + [90]*3
        cog = self.atoms.center_of_geometry(wrap=False)
        box_center = box / 2
        self.atoms.translate(box_center - cog)

    def _write(self, selection, name):
        atoms = self.polymer.select_atoms(selection)
        atoms.write(name)

    def save(self, dummies="X*", fname="polymer", selectionString = None, gmx = False):
        """
        Save polymer as a PDB or GROMACS file with dummy atoms excluded.
        Optionally, select a subset of the polymer to save.

        :param dummies: names of all the dummy atoms, for use in the 
                    selection string (e.g. 'CN CMA CP CQ' to exclude these 
                    four dummy atom types), defaults to "X*"
        :type dummies: str, optional
        :param fname: name of the output file, defaults to "polymer"
        :type fname: str, optional
        :param selectionString: MDAnalysis atom selection string, for more 
                    information on supported atom selection formats see 
                    `MDAnalysis Atom selection language <https://userguide.mdanalysis.org/stable/selections.html>`, 
                    defaults to None
        :type selectionString: _type_, optional
        :param gmx: save output as GROMACS if True, else save as default PDB, 
                defaults to False
        :type gmx: bool, optional
        """
        if selectionString:
            if gmx:
                self.select_atoms(f"{selectionString} and not name {dummies}").atoms._write(f"{fname}.gro")
            else:
                self.select_atoms(f"{selectionString} and not name {dummies}").atoms._write(f"{fname}.pdb")
        else:
            if gmx:
                self._write(f"not name {dummies}", f"{fname}.gro")
            else:
                self._write(f"not name {dummies}", f"{fname}.pdb")

    def crudesave(self,fname="polymer_crude"):
        """
        Save polymer as a PDB file with dummy atoms included.
        Useful for debugging the geometry transforms.

        :param fname: name of the output file, defaults to "polymer_crude"
        :type fname: str, optional
        """
        self.atoms._write(f"{fname}.pdb")

    def select_atoms(self, selection) -> mda.AtomGroup:
        """
        Selection method that selects atoms from the polymer's MDAnalysis
        Universe by leveraging the MDAnalysis atom selection.

        :param selection: MDAnalysis atom selection string, for more 
                    information on supported atom selection formats see 
                    `MDAnalysis Atom selection language <https://userguide.mdanalysis.org/stable/selections.html>`
        :type selection: str
        :return: An MDAnalysis AtomGroup containing only the selected atoms
        :rtype: mda.AtomGroup
        """
        return self.polymer.select_atoms(selection)
