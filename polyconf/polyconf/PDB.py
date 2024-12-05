#!/usr/bin/env python 

#TODO we need to change the name of this if the default is not to save as a pdb
    
class PDB: 
    """
    docstring go brr
    """
    def __init__(self, polymer) -> None:
        """
        Save a copy of the Polymer.polymer MDAnalysis Universe and its atoms
        for cleanup and PDB export.

        Args:
            polymer (MDAnalysis Universe): polymer constructed with the Polymer
                    class, supplied from a Polymer.polymer attribute
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

    def save(self, dummyAtoms="X*", fname="polymer", selectionString = None, gmx = False):
        """
        Save polymer as a GROMACS .gro file with dummy atoms excluded.
        Optionally, select a subset of the polymer to save.

        Args:
            dummyAtoms (str): names of all the dummy atoms, for use in the 
                    selection string (e.g. 'CN CMA CP CQ' to exclude these 
                    four dummy atom types)
            fname (str): name of the output file, default is 'polymer'
            selectionString (str): MDAnalysis atom selection string, for more 
                    information on supported atom selection formats see 
                    [MDAnalysis Atom selection language](https://userguide.mdanalysis.org/stable/selections.html)
            gmx (bool): save output as GROMACS if True, else save as default PDB
        """
        if selectionString:
            if gmx:
                self.select_atoms(f"{selectionString} and not name {dummyAtoms}").atoms._write(f"{fname}.gro")
            else:
                self.select_atoms(f"{selectionString} and not name {dummyAtoms}").atoms._write(f"{fname}.pdb")
        else:
            if gmx:
                self._write(f"not name {dummyAtoms}", f"{fname}.gro")
            else:
                self._write(f"not name {dummyAtoms}", f"{fname}.pdb")

    def crudesave(self,fname="polymer_crude"):
        """
        Save polymer as a PDB file with dummy atoms included.
        Useful for debugging the geometry transforms.

        Args:
            fname (str): name of the output file, default is 'polymer_crude'
        """
        self.atoms._write(f"{fname}.pdb")

    def select_atoms(self, selection):
        """
        Selection method that selects atoms from the polymer's MDAnalysis
        Universe by leveraging the MDAnalysis atom selection.

        Args:
            selection (str): MDAnalysis atom selection string, for more 
                    information on supported atom selection formats see 
                    [MDAnalysis Atom selection language](https://userguide.mdanalysis.org/stable/selections.html)
        """
        return self.polymer.select_atoms(selection)
