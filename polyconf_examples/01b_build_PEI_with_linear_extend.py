#!/usr/bin/env python

from polyconf.Monomer import Monomer
from polyconf.Polymer import Polymer
from polyconf.PDB import PDB


###############################################################################
#                                                                             #
#  Polymer two: Building a linear polymer with linear extend()                #
#                                                                             #
###############################################################################

# this polymer will also be 128 monomers long, and we will use the same input files as before
# however, this time we will extend the polymer by along a vector [1,0,0].
# this will result in a pseudolinear conformation, which should greatly reduce the number of overlapping atoms

polymer2=Polymer(Monomer('PEI_start.pdb')) # initialise the polymer the same way

adds=126 
aligner=[1,0,0] # linearization vector

for _ in range (0,adds):

    polymer2.extend(
            Monomer('PEI_monomer.pdb'),
            n=polymer2.maxresid(), 
            nn=polymer2.newresid(),
            names=dict(P='C1',Q='CX',R='NX',S='N1', V1='C1',V2='N1'), # we are aligning each monomer so the vector A1 A2 lies along aligner
            linearise=True, # linearize incoming monomer so monomer vector lies along a specified vector, named ortho
            ortho=aligner,
            joins=[('N1','C1')], 
            ) 

# now we add the final monomer
# the process is the same, but we are using a different pdb file

polymer2.extend(
            Monomer('PEI_end.pdb'),
            n=polymer2.maxresid(), 
            nn=polymer2.newresid(),
            names=dict(P='C1',Q='CX',R='NX',S='N1', V1='C1',V2='N1'), # we are aligning each monomer so the vector A1 A2 lies along aligner
            linearise=True, # linearize incoming monomer so monomer vector lies along along ortho
            ortho=aligner,
            joins=[('N1','C1')], 
            ) 


# you may have noticed some differences

# first, we did not define n and new_n separately
# we just called maxresid() and newresid() during extend()

# we also didn't use renamer to flag dummy atoms
# this is because our monomer coordinate files were designed so every monomer has the same atoms names, 
# so that our dummy atoms would be be every atom named NX or CX 

# when we save our polymer, we can explicitly name our dummy atoms

Saver = PDB(polymer2)
Saver.cleanup()
Saver.save(dummies='CX NX',fname='PEI_linear_02_simple-linear-extend') # save, excluding dummy atoms

# When you examine polymer two , you can see that the resulting structure is extended in a long linear chain.  
# The linearization process has introduced some unrealistic bond angles, the angle [ N1_n , C1_n+1 , C2_n+1] is completely linear.  
# Additionally, the structure is still highly ordered.

# These pseudolinear structures can be useful for inspecting your polymer to ensure you have built it correctly. 
# They may also be useful starting points for generating valid conformations with shuffle() and dihedral_solver()
# If your simulation parameters are correct, the unrealistic angles will be corrected during energy minimzation and equilibration. 

###############################################################################
#                                                                             #
#  Polymer three: Building a linear polymer with alternating extend() steps   #
#                                                                             #
###############################################################################

# this polymer will also be 128 monomers long, and we will use the same input files as before
# however, this time we will align to alternating vectors
# even monomers will be aligned to the vector [1,0,0]
# odd monomers will be aligned to the vector [0,1,0] 

# they will continue to switch between the two until they reach the end

polymer3=Polymer(Monomer('PEI_start.pdb')) # initialise
adds=126 # lay 24 additional monomers
alternator=True # this is a boolean we will use to alternate between our two alignment vectors

vector1=[1,0,0] # linearization vector
vector2=[0,1,0] # linearization vector

for i in range (0,adds):
    # first, choose which alignment vector to use

    if alternator:
        aligner=vector1
    else: 
        aligner=vector2

    alternator=not alternator #flip the alternator

    polymer3.extend(
            Monomer('PEI_monomer.pdb'),
            n=polymer3.maxresid(), 
            nn=polymer3.newresid(),
            names=dict(P='C1',Q='CX',R='NX',S='N1', V1='C1',V2='N1'), 
            linearise=True,
            ortho=aligner,
            joins=[('N1','C1')], 
            ) 

# then add the final monomer

    if alternator:
        aligner=vector1
    else: 
        aligner=vector2

polymer3.extend(
            Monomer('PEI_end.pdb'),
            n=polymer3.maxresid(), 
            nn=polymer3.newresid(),
            names=dict(P='C1',Q='CX',R='NX',S='N1', V1='C1',V2='N1'), 
            linearise=True,
            ortho=aligner,
            joins=[('N1','C1')], 
            ) 

# save this polymer

Saver = PDB(polymer3)
Saver.cleanup() # center in box
Saver.save(dummies='CX NX',fname='PEI_linear_03_alternating-linear-extend') # save, excluding dummy atoms

# When you examine polymer three , you can see that the resulting structure is again extended in a long linear chain.  
# The conformation of this chain is different, and the bond angles are less unrealistic.  If your simulation parameters are correct, energy minimzation would ajust the bond angles
# however there are still some problems.
# This conformation has ammonium groups positioned close to each other, despite the repulsive electrostatic interactions betweeen them.
# Additionally, the structure is still highly ordered.

# For our next polymer, we will address the problem of highly ordered starting structures.
