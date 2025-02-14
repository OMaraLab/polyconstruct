from polyconf.Polymer import Polymer
from polyconf.Monomer import Monomer
from polyconf.PDB import PDB

polymer = Polymer(Monomer("monomers/UNK_460A12_PolyConf.pdb"))

for i in range(29):
    polymer.extend(Monomer("monomers/UNK_460A12_PolyConf.pdb"), n=polymer.maxresid(), nn=polymer.newresid(),
            names=dict(P1='C04',P2='C02',Q1='C03',Q2='C01',R='C02',S='C03'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i
            joins=[('C02','C03')],
            )
    polymer.renamer(polymer.maxresid()-1, "C04", nameout='X')
    polymer.renamer(polymer.maxresid()-1, "H0K", nameout='X')
    polymer.renamer(polymer.maxresid()-1, "H0M", nameout='X')
    polymer.renamer(polymer.maxresid()-1, "C05", nameout='X')
    polymer.renamer(polymer.maxresid()-1, "H0N", nameout='X')
    polymer.renamer(polymer.maxresid()-1, "H0O", nameout='X')
    polymer.renamer(polymer.maxresid()-1, "H0P", nameout='X')

    polymer.renamer(polymer.maxresid(), "C01", nameout='X')
    polymer.renamer(polymer.maxresid(), "H0F", nameout='X')
    polymer.renamer(polymer.maxresid(), "H0G", nameout='X')
    polymer.renamer(polymer.maxresid(), "C00", nameout='X')
    polymer.renamer(polymer.maxresid(), "H0C", nameout='X')
    polymer.renamer(polymer.maxresid(), "H0D", nameout='X')
    polymer.renamer(polymer.maxresid(), "H0E", nameout='X')

Saver = PDB(polymer)
Saver.cleanup() # center the Polymer in the PBC box
Saver.save(dummyAtoms='X*',fname='polymers/30mer')
