#!/usr/bin/env python

from polyconf.Monomer import Monomer
from polyconf.Polymer import Polymer
from polyconf.PDB import PDB

import random

random.seed(12) # by setting the random seed we can ensure the random conformations are replicable.  

###############################################################################
#                                                                             #
#  Polymer six: Generating an ensemble of conformations                       #
#                                                                             #
###############################################################################

# In this example, we will follow the process used in polymer 5 to generate a conformation,
# but we will repeat it five times to generate an ensemble of starting conformations.

# Using different starting conformations reduces sampling bias arising from the starting conformation.
# This is very useful when you wish to study polymer dynamics

# First, make a PEI polymer through the same process as polymer 1

polymer6=Polymer(Monomer('PEI_start.pdb')) 

adds=126 

for i in range (0,adds):
    polymer6.extend( # extend with one monomer, aligned along this step's linearization vector
            Monomer('PEI_monomer.pdb'), # extend with this monomer
            n=polymer6.maxresid(), # we will allways add onto the existing monomer with the highest resid
            nn=polymer6.newresid(), # the incoming monomer needs a new resid
            names=dict(P='C1',Q='CX',R='NX',S='N1'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
            joins=[('N1','C1')],# new connection between N1_i and C1_i+1 
            ) 
polymer6.extend(
            Monomer('PEI_end.pdb'),
            n=polymer6.maxresid(), 
            nn=polymer6.newresid(),
            names=dict(P='C1',Q='CX',R='NX',S='N1'),
            joins=[('N1','C1')],
            )


# Generate lists of dihedrals as for polymer 5

NC_dihedrals=polymer6.gen_pairlist(J='N1',K='C1',first_resid=1,last_resid=127,mult=3,same_res=False) 
alkane_dihedrals=polymer6.gen_pairlist(J='C1',K='C2',first_resid=1,last_resid=128,mult=3)

# Then, we will generate five starting conformations

# For each conformation:
#   start by making a copy of our initial conformation
#   randomise the conformation by shuffling the NC dihedrals
#   solve the conformation by rotating the alkane dihedrals

for conf in range(1,6): 

    print(f'generating conformation {conf}')
    newconf=polymer6.copy() # make a duplicate of the original polymer

    newconf.shuffler(NC_dihedrals)  # as before, randomise the NC dihedrals to generate a random conformation 
    newconf.dihedral_solver(alkane_dihedrals,dummies='CX NX',cutoff=0.9) # solve by rotating the alkane dihedrals

    Saver = PDB(newconf)
    Saver.cleanup()
    Saver.save(dummies='CX NX',fname=f'PEI_linear_06_conformation_{conf}') 

# This should produce an ensemble of five starting conformations
