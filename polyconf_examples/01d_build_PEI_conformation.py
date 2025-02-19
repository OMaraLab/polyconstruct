#!/usr/bin/env python

from polyconf.Monomer import Monomer
from polyconf.Polymer import Polymer
from polyconf.PDB import PDB


###############################################################################
#                                                                             #
#  Polymer five: Generating a randomised starting conformation shuffler()     #
#                                                                             #
###############################################################################

# in this example, we will create a PEI polymer as in example 01a, then use shuffler() to generate a random conformation
# we will then use dihedral_solver() to convert the random conformation into something reasonable

#first, make a copy of polymer one

polymer5=Polymer(Monomer('PEI_start.pdb')) 

# then, we will add the middle monomers
# we will do this by adding one middle monomer, and repeating 126 times 

adds=126 

for i in range (0,adds):
    polymer5.extend( # extend with one monomer, aligned along this step's linearization vector
            Monomer('PEI_monomer.pdb'), # extend with this monomer
            n=polymer5.maxresid(), # we will allways add onto the existing monomer with the highest resid
            nn=polymer5.newresid(), # the incoming monomer needs a new resid
            names=dict(P='C1',Q='CX',R='NX',S='N1'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
            joins=[('N1','C1')],# new connection between N1_i and C1_i+1 
            ) 

# now we add the final monomer
# the process is the same, but we are using a different pdb file

polymer5.extend(
            Monomer('PEI_end.pdb'),
            n=polymer5.maxresid(), 
            nn=polymer5.newresid(),
            names=dict(P='C1',Q='CX',R='NX',S='N1'),
            joins=[('N1','C1')],
            )



import random

random.seed(8) # by setting the random seed we can ensure the random conformations are replicable.  

# we are going to generate a conformation in a two step process
# first, we will are going to randomly rotate the dihedrals centered on the bond between atoms N1_n and C1_n+1
# then, we are going to solve the polymer by rotating each dihedral centered on a C1-C2 bond

# in principle, you can shuffle as many dihedrals as you like
# however, the process of randomly rotating dihedrals does not consider interatomic distances
# in my experience, randomising more dihedrals results in very tight compacted conformations with may overlapping atoms
# it is best to choose only one ot two bonds to shuffle
# if you are not having success shuffling one bond, it may be useful to try shuffling a different bond instead
# or shuffling a second time to find a more amenable starting conformation

NC_dihedrals=polymer5.gen_pairlist(J='N1',K='C1',first_resid=1,last_resid=127,mult=3,same_res=False) 
# these dihedrals are centered on the bond that links two monomers, so we need to set same_res to False

alkane_dihedrals=polymer5.gen_pairlist(J='C1',K='C2',first_resid=1,last_resid=128,mult=3)

polymer5.shuffler(NC_dihedrals,)  # this generates a random conformation 

# lets look at the intermediate conformation we got by random shuffling

Saver = PDB(polymer5)
Saver.cleanup() # center in box
Saver.save(dummies='CX NX',fname='PEI_linear_05_shuffled') # save, excluding dummy atoms

# now solve the conformation
polymer5.dihedral_solver(alkane_dihedrals,dummies='CX NX',cutoff=0.9) # this converts the shuffled conformation into one without overlapping atoms, by rotating the C1-C2 torsions

Saver = PDB(polymer5)
Saver.cleanup() # center in box
Saver.save(dummies='CX NX',fname='PEI_linear_05_shuffled_then_solved') # save, excluding dummy atoms

# Polymer five has two output files, shuffled, and shuffled_then_solved. 
# When you examine polymer five, you can see that the two conformations are extended random coils, and is quite different to the conformations of polymers one to four.

# The shuffled conformation has some steric clashes.  For example, monomers 47 and 54 are overlapping.
# The shuffled_then_solved conformation has been adjusted to remove these clashes.

