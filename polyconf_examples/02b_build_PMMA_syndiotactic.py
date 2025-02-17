#!/usr/bin/env python

from polyconf.Monomer import Monomer
from polyconf.Polymer import Polymer
from polyconf.PDB import PDB
import random

random.seed(1) 


###############################################################################
#                                                                             #
#  Polymer two: Building a syndiotactic PMMA polymer                          #
#                                                                             #
###############################################################################

print ('creating a syndiotactic PMMA polymer') 

dummies="CMA CN CP CQ"

Monomer_D='MMAD_bonds.pdb'
Monomer_L='MMAL_bonds.pdb'
alternator=True

PMMA_syndiotactic=Polymer(Monomer(Monomer_D)) # initialise

adds=49 # we will lay 49 additional monomers

for i in range (0,adds):
        if alternator:
            monomer=Monomer_L
        else:
            monomer=Monomer_D
        alternator=not alternator
        PMMA_syndiotactic.extend( # extend with one monomer, aligned along this step's linearization vector
            Monomer(monomer), # extend with this monomer
            n=PMMA_syndiotactic.maxresid(), # extend existing residue i
            nn=PMMA_syndiotactic.newresid(), # incoming monomer will have resid i+1
            names=dict(Q='CA',P='CMA',S='C',R='CN'), # C1_i+1 fit to CX_i
            joins=[('C','CA')],# new connection between N1_i and C1_i+1 
            ) 

Saver = PDB(PMMA_syndiotactic)
Saver.cleanup() # center in box
Saver.save(dummies=dummies,fname='PMMA_syndiotactic_vanilla-extend') # save, excluding dummies atoms

alkanes=PMMA_syndiotactic.gen_pairlist(J='C',K='CA',first_resid=1,same_res=False,last_resid=49,mult=3)
sidechains=PMMA_syndiotactic.gen_pairlist(J='CA',K='CB',first_resid=1,same_res=True,last_resid=50,mult=6)

PMMA_syndiotactic.dihedral_solver(alkanes,dummies=dummies,cutoff=1.1) # this converts the shuffled conformation into one without overlapping atoms
PMMA_syndiotactic.dihedral_solver(sidechains,dummies=dummies,cutoff=1.1) # this converts the shuffled conformation into one without overlapping atoms

Saver = PDB(PMMA_syndiotactic)
Saver.cleanup() # center in box
Saver.save(dummies=dummies,fname='PMMA_syndiotactic_vanilla-solved') # save, excluding dummies atoms


PMMA_syndiotactic.shuffler(alkanes)
PMMA_syndiotactic.dihedral_solver(alkanes,dummies=dummies,cutoff=1.1) # this converts the shuffled conformation into one without overlapping atoms
PMMA_syndiotactic.shuffler(sidechains)
PMMA_syndiotactic.dihedral_solver(sidechains,dummies=dummies,cutoff=1.1) # this converts the shuffled conformation into one without overlapping atoms

Saver = PDB(PMMA_syndiotactic)
Saver.cleanup() # center in box
Saver.save(dummies=dummies,fname='PMMA_syndiotactic_shuffled_and_solved') # save, excluding dummies atoms


# shuffler() will generate a random conformation by shuffling dihedrals.  
# Some conformations generated with shuffler() cannot be easilly solved.  In this example, we have chosen a random seed that will work reliably.
# Controlling the random seed can also be useful for replicability in your scripts.
# After you have completed the tutorial, try varying the random seed, and comparing the resulting conformations.
