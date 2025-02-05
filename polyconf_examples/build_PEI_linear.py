#!/usr/bin/env python

from polyconf.Monomer import Monomer
from polyconf.Polymer import Polymer
from polyconf.PDB import PDB


###############################################################################
#                                                                             #
#  Polymer one: Building a polymer with  extend()                             #
#                                                                             #
###############################################################################

polymer1=Polymer(Monomer('PEI_start.pdb')) # initialise
imax=127 # we will lay 127 additional monomers

for i in range (0,imax):
    if not i==imax: # do this for everything except the first one
        polymer1.extend( # extend with one monomer, aligned along this step's linearization vector
            Monomer('PEI_monomer.pdb'), # extend with this monomer
            n=polymer1.maxresid(), # extend existing residue i
            nn=polymer1.newresid(), # incoming monomer will have resid i+1
            names=dict(P1='CX',P2='C1',Q1='N1',Q2='NX'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
            joins=[('N1','C1')],# new connection between N1_i and C1_i+1 
            ) 
    else:
        # add final monomer
        polymer1.extend(Monomer('PEI_end.pdb'),i,i+1,
            names=dict(P1='CX',P2='C1',Q1='N1',Q2='NX'),
            joins=[('N1','C1')],)

Saver = PDB(polymer1)
Saver.cleanup() # center in box
Saver.save(dummyAtoms='CX NX',fname='polymer_01_vanilla-extend') # save, excluding dummy atoms

# When you examine polymer one , you can see that the resulting strucure is a tightly coiled helix.  
# This structure is highly ordered, and the turns of the helix are very close.



###############################################################################
#                                                                             #
#  Polymer two: Building a linear polymer with linear extend()                #
#                                                                             #
###############################################################################


polymer2=Polymer(Monomer('PEI_start.pdb')) # initialise
imax=127 # lay 127 additional monomers
ortho=[1,0,0] # linearization vector

for i in range (1,imax+1):
    if not i==imax: # do this for everything except the first one
        polymer2.extend(
            Monomer('PEI_monomer.pdb'), # extend with this monomer
            i, # extend existing residue i
            i+1, # incoming monomer will have resid i+1
            names=dict(P1='CX',P2='C1',Q1='N1',Q2='NX', # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
            R='C1',S='N1'), # monomer vector is vector C1_i+1 N1_i+1,
            linearise=True, # linearize incoming monomer so monomer vector lies along along ortho
            ortho=ortho,# specify linearization vector
            joins=[('N1','C1')],# new connection between N1_i and C1_i+1 
            ) 
    else:
        # add final monomer
        polymer2.extend(Monomer('PEI_end.pdb'),i,i+1,
            names=dict(P1='CX',P2='C1',Q1='N1',Q2='NX'),
            joins=[('N1','C1')],)
#    polymer.genconf([('C2','N1')],length=10, cutoff = 0.5)
    'C1' 'NX' 'CX' 'N1'

Saver = PDB(polymer2)
Saver.cleanup() # center in box
Saver.save(dummyAtoms='CX NX',fname='polymer_02_simple-linear-extend') # save, excluding dummy atoms

# When you examine polymer two , you can see that the resulting structure is extended in a long linear chain.  
# The linearization process has introduced some unrealistic bond angles, the angle [ N1_n , C1_n+1 , C2_n+1] is linear.  
# Additionally, the structure is still highly ordered.

# These pseudolinear structures can be useful for inspecting your polymer to ensure you have built it correctly. 
# They may also be useful starting points for generating valid conformations with shuffle() and dihedral_solver()
# If your simualtion pramaters are correct, the unrealistic angles will be corrected during energy minimzation and equilibration. 


###############################################################################
#                                                                             #
#  Polymer three: Building a linear polymer with alternating extend() steps   #
#                                                                             #
###############################################################################

# In this example, we will use a linear extend, but we will alternate between two vectors
# the first linear extend will use vector A
# the second linear extend will use vector B
# they will continue to switch between the two until they reach the end

polymer3=Polymer(Monomer('PEI_start.pdb')) # initialise
imax=127 # lay 24 additional monomers
alternator=True # lay monomers along alternating vectors to get a clean linear conformation, because shuffle isn't working



for i in range (1,imax+1):
    if not i==imax: # do this for everything except the first one
        if alternator:
            ortho=[1,0,0] # linearization vector A
        else:
            ortho=[0,1,0] # linearization vector B
        alternator=not alternator # flip the alternator

        polymer3.extend(Monomer( # extend with one monomer, aligned along this step's linearization vector
            'PEI_monomer.pdb'), # extend with this monomer
            i, # extend existing residue i
            i+1, # incoming monomer will have resid i+1
            names=dict(P1='CX',P2='C1',Q1='N1',Q2='NX', # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
            R='C1',S='N1'), # monomer vector is vector C1_i+1 N1_i+1,
            linearise=True, # linearize incoming monomer so monomer vector lies along along ortho
            ortho=ortho,# specify linearization vector
            joins=[('N1','C1')],# new connection between N1_i and C1_i+1 
            ) 
    else:
        # add final monomer
        polymer3.extend(Monomer('PEI_end.pdb'),i,i+1,
            names=dict(P1='CX',P2='C1',Q1='N1',Q2='NX'),
            joins=[('N1','C1')],)
#    polymer.genconf([('C2','N1')],length=10, cutoff = 0.5)
    'C1' 'NX' 'CX' 'N1'

Saver = PDB(polymer3)
Saver.cleanup() # center in box
Saver.save(dummyAtoms='CX NX',fname='polymer_03_alternating-linear-extend') # save, excluding dummy atoms

# When you examine polymer three , you can see that the resulting structure is again extended in a long linear chain.  
# The conformation of this chain is different, and the bond angles are more realistic.
# This conformation has ammonium groupos positioned close to each other, and may be high energy.
# Additionally, the structure is still highly ordered.


###############################################################################
#                                                                             #
#  Polymer four: Building a reasonable conformation with dihedral_solver()    #
#                                                                             #
###############################################################################

# in this example, we will make a copy of polymer 1, then use dihedral_solver() to generate a more reasonable starting conformation

polymer4=Polymer(Monomer('PEI_start.pdb')) # initialise
imax=127 # lay 127 additional monomers

# TODO replace with deepcopy

for i in range (1,imax+1):
    if not i==imax: # do this for everything except the first one
        polymer4.extend( # extend with one monomer, aligned along this step's linearization vector
            Monomer('PEI_monomer.pdb'), # extend with this monomer
            i, # extend existing residue i
            i+1, # incoming monomer will have resid i+1
            names=dict(P1='CX',P2='C1',Q1='N1',Q2='NX'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
            joins=[('N1','C1')],# new connection between N1_i and C1_i+1 
            ) 
    else:
        # add final monomer
        polymer4.extend(Monomer('PEI_end.pdb'),i,i+1,
            names=dict(P1='CX',P2='C1',Q1='N1',Q2='NX'),
            joins=[('N1','C1')],)


CN=polymer4.gen_pairlist(a1='C2',a2='N1',first_resid=1,last_resid=128,mult=3)


polymer4.dihedral_solver(CN,dummy='CX,NX',cutoff=0.7)
# this will iterate through every bond in C2-N1 bond, rotating them so that atoms one either side of the bond are always more than 0.7 A apart.

Saver = PDB(polymer4)
Saver.cleanup() # center in box
Saver.save(dummyAtoms='CX NX',fname='polymer_04_solved') # save, excluding dummy atoms

# When you examine polymer four, you can see that the resulting structure is much more reasonable than polymer 1.  
# The polymer is arranged in a twisting helix.  The loops are much less tightly packed, and the distances between adjacent loops are much more reasonable.
# Furthermore, because we did not use a linear extend, the bonds and angles in this structure are consistent with those in the original monomers
# However, the algorithmic solver has still resulted in a highly ordered structure.  
# This polymer may be reasonable for illustrative purposes, but we may wish to generate an ensemble of less ordered structure for production simulations.




###############################################################################
#                                                                             #
#  Polymer five: Generating a randomised starting conformation shuffler()     #
#                                                                             #
###############################################################################

# in this example, we will make a copy of polymer 1, then use shuffler() to generate a random conformation
# we will then use dihedral_solver() to convert the random conformation into something reasonable

#first, make a copy of polymer one

polymer5=Polymer(Monomer('PEI_start.pdb')) # initialise
imax=127 # lay 127 additional monomers


for i in range (1,imax+1):
    if not i==imax: # do this for everything except the first one
        polymer5.extend( # extend with one monomer, aligned along this step's linearization vector
            Monomer('PEI_monomer.pdb'), # extend with this monomer
            i, # extend existing residue i
            i+1, # incoming monomer will have resid i+1
            names=dict(P1='CX',P2='C1',Q1='N1',Q2='NX'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
            joins=[('N1','C1')],# new connection between N1_i and C1_i+1 
            ) 
    else:
        # add final monomer
        polymer5.extend(Monomer('PEI_end.pdb'),i,i+1,
            names=dict(P1='CX',P2='C1',Q1='N1',Q2='NX'),
            joins=[('N1','C1')],)


import random

random.seed(8) # by setting the random seed we can ensure the random conformations are replicable

alkanes=polymer5.gen_pairlist(a1='C1',a2='C2',first_resid=1,last_resid=128,mult=3)
NC=polymer5.gen_pairlist(a1='N1',a2='C1',first_resid=1,last_resid=127,mult=6,same_res=False)

polymer5.shuffler(NC,)  # this generates a random conformation 

Saver = PDB(polymer5)
Saver.cleanup() # center in box
Saver.save(dummyAtoms='CX NX',fname='polymer_05_shuffled') # save, excluding dummy atoms


polymer5.dihedral_solver(alkanes,dummy='CX NX',cutoff=1) # this converts the shuffled conformation into one without overlapping atoms


Saver = PDB(polymer5)
Saver.cleanup() # center in box
Saver.save(dummyAtoms='CX NX',fname='polymer_05_shuffled_then_solved') # save, excluding dummy atoms

# Polymer five has two output files, shuffled, and shuffled_then_solved. 
# When you examine polymer five, you can see that the two conformations are extended random coils, and is quite different to the conformations of polymers one to four.

# The shuffled conformation has some steric clashes.  For example, monomers 47 and 54 are overlapping.
# The shuffled_then_solved conformation has been adjusted to remove these clashes.