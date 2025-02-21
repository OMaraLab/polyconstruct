#!/usr/bin/env python

from polyconf.Monomer import Monomer
from polyconf.Polymer import Polymer
from polyconf.PDB import PDB
import random
random.seed(1) 

###############################################################################
#                                                                             #
#  Polymer one: Building an isotactic PMMA polymer                            #
#                                                                             #
###############################################################################

print ('creating an isotactic PMMA polymer') 

dummies="CMA CN CP CQ"

PMMA_isotactic=Polymer(Monomer('MMAD_bonds.pdb')) # initialise

adds=49 # we will lay 49 additional monomers

for i in range (0,adds):
        PMMA_isotactic.extend( # extend with one monomer, aligned along this step's linearization vector
            Monomer('MMAD_bonds.pdb'), # extend with this monomer
            n=PMMA_isotactic.maxresid(), # extend existing residue i
            nn=PMMA_isotactic.newresid(), # incoming monomer will have resid i+1
            names=dict(P='CMA',Q='CA',R='CN',S='C'), # CA_i+1 fit to CMA_i, then rotate so C_i+1 fit to CN_i 
            joins=[('C','CA')],# new connection between C_i and CA_i+1 
            ) 

Saver = PDB(PMMA_isotactic)
Saver.cleanup() # center in box
Saver.save(dummies=dummies,fname='PMMA_isotactic_vanilla-extend') # save, excluding dummy atoms

# next, generate conformations by adjusting the dihedrals

alkanes=PMMA_isotactic.gen_pairlist(J='C',K='CA',first_resid=1,same_res=False,last_resid=49,mult=3)
sidechains=PMMA_isotactic.gen_pairlist(J='CA',K='CB',first_resid=1,same_res=True,last_resid=50,mult=6)

PMMA_isotactic.dihedral_solver(alkanes,dummies=dummies,cutoff=0.8) # this converts the shuffled conformation into one without overlapping atoms by rotating the alkanes

Saver = PDB(PMMA_isotactic)
Saver.cleanup() # center in box
Saver.save(dummies=dummies,fname='PMMA_isotactic_vanilla-solved') 

PMMA_isotactic.shuffler(alkanes) # randomise the rotation of C-CA torsions
PMMA_isotactic.dihedral_solver(alkanes,dummies=dummies,cutoff=1.1) # this converts the shuffled conformation into one without overlapping atoms

PMMA_isotactic.shuffler(sidechains) # randomise the position of the sidechains
PMMA_isotactic.dihedral_solver(sidechains,dummies=dummies,cutoff=0.8) 

Saver = PDB(PMMA_isotactic)
Saver.cleanup() # center in box
Saver.save(dummies=dummies,fname='PMMA_isotactic_shuffled_and_solved') 


