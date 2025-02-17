#!/usr/bin/env python

from polyconf.Monomer import Monomer
from polyconf.Polymer import Polymer
from polyconf.PDB import PDB


###############################################################################
#                                                                             #
#  Polymer one: Building a polymer with  extend()                             #
#                                                                             #
###############################################################################

# this polymer will be 128 monomers long.
# there will be an initial monomer ( PEI_start.pdb ), 126 middle monomers (PEI_monomer.pdb) , and a final monomer ('PEI_end.pdb')

# first, we initialise a new polymer with the starting monomer

polymer1=Polymer(Monomer('PEI_start.pdb')) 

# then, we will add the middle monomers
# we will do this by adding one middle monomer, and repeating 126 times 

adds=126 

for i in range (0,adds):
    polymer1.extend( # extend with one monomer, aligned along this step's linearization vector
            Monomer('PEI_monomer.pdb'), # extend with this monomer
            n=polymer1.maxresid(), # we will allways add onto the existing monomer with the highest resid
            nn=polymer1.newresid(), # the incoming monomer needs a new resid
            names=dict(Q='CX',P='C1',S='N1',R='NX'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
            joins=[('N1','C1')],# new connection between N1_i and C1_i+1 
            ) 


# now we add the final monomer
# the process is the same, but we are using a different pdb file

polymer1.extend(
            Monomer('PEI_end.pdb'),
            n=polymer1.maxresid(), 
            nn=polymer1.newresid(),
            names=dict(Q='CX',P='C1',S='N1',R='NX'),
            joins=[('N1','C1')],
            )

# Now we can save our polymer as a pdb

# we start by making a PDB object, which will save our file

Saver = PDB(polymer1)

Saver.cleanup() # this puts our polymer in the center of the box

# saver can automatically remove dummy atoms
# here we are saying that every atom with the name CX or NX is a dummy atom
# this works because we set up our files so those dummy atoms had consistent names in every single monomer

Saver.save(dummyAtoms='CX NX',fname='polymer_01_vanilla-extend') 

# When you examine polymer one , you can see that the resulting strucure is a tightly coiled helix.  
# This structure is highly ordered, and the turns of the helix are very close.
# It's not a very likely conformation

# Lets use some other methods to make different conformations

###############################################################################
#                                                                             #
#  Polymer two: Building a linear polymer with linear extend()                #
#                                                                             #
###############################################################################

# this polymer will also be 128 monomers long, and we will use the same input files as before
# however, this time we will extend the polymer by along a vector [1,0,0].

polymer2=Polymer(Monomer('PEI_start.pdb')) # initialise the polymer the same way

adds=126 
aligner=[1,0,0] # linearization vector

for i in range (0,adds):

    polymer2.extend(
            Monomer('PEI_monomer.pdb'),
            n=polymer2.maxresid(), 
            nn=polymer2.newresid(),
            names=dict(Q='CX',P='C1',S='N1',R='NX', V1='C1',V2='N1'), # we are aligning each monomer so the vector A1 A2 lies along aligner
            linearise=True, # linearize incoming monomer so monomer vector lies along along ortho
            ortho=aligner,
            joins=[('N1','C1')], 
            ) 

# now we add the final monomer
# the process is the same, but we are using a different pdb file
polymer2.extend(
            Monomer('PEI_end.pdb'),
            n=polymer2.maxresid(), 
            nn=polymer2.newresid(),
            names=dict(Q='CX',P='C1',S='N1',R='NX', V1='C1',V2='N1'), # we are aligning each monomer so the vector A1 A2 lies along aligner
            linearise=True, # linearize incoming monomer so monomer vector lies along along ortho
            ortho=aligner,
            joins=[('N1','C1')], 
            ) 


# save the polymer as before

Saver = PDB(polymer2)
Saver.cleanup() # center in box
Saver.save(dummyAtoms='CX NX',fname='polymer_02_simple-linear-extend') # save, excluding dummy atoms

# When you examine polymer two , you can see that the resulting structure is extended in a long linear chain.  
# The linearization process has introduced some unrealistic bond angles, the angle [ N1_n , C1_n+1 , C2_n+1] is completely linear.  
# Additionally, the structure is still highly ordered.

# These pseudolinear structures can be useful for inspecting your polymer to ensure you have built it correctly. 
# They may also be useful starting points for generating valid conformations with shuffle() and dihedral_solver()
# If your simualtion pramaters are correct, the unrealistic angles will be corrected during energy minimzation and equilibration. 


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
            names=dict(Q='CX',P='C1',S='N1',R='NX', V1='C1',V2='N1'), 
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
            names=dict(Q='CX',P='C1',S='N1',R='NX', V1='C1',V2='N1'), 
            linearise=True,
            ortho=aligner,
            joins=[('N1','C1')], 
            ) 

# save this polymer

Saver = PDB(polymer3)
Saver.cleanup() # center in box
Saver.save(dummyAtoms='CX NX',fname='polymer_03_alternating-linear-extend') # save, excluding dummy atoms

# When you examine polymer three , you can see that the resulting structure is again extended in a long linear chain.  
# The conformation of this chain is different, and the bond angles are less unrealistic.  Hopefully, energy minimzation would correct this
# however there are still some problems.
# This conformation has ammonium groupos positioned close to each other, and may be high energy.
# Additionally, the structure is still highly ordered.

# for our next polymer, we will address the problem of highly ordered starting structures.

###############################################################################
#                                                                             #
#  Polymer four: Building a reasonable conformation with dihedral_solver()    #
#                                                                             #
###############################################################################

# in this example, we will make a copy of polymer 1, then use dihedral_solver() to generate a more reasonable starting conformation

polymer4=polymer1.copy() # make a copy of our first polymer, which was a tightly packed helix 

# we need to choose a set of deihdrals to shuffle
# I am going to choose every dihedral centered on the bond between atoms C2 and N1

# we use gen pairlist to generate a list of all of these dihedrals
# note that we only specify the two atoms in the center of the dihedral

CN_dihedrals=polymer4.gen_pairlist(J='C2',K='N1',first_resid=1,last_resid=128,mult=3) 

# next, we use dihedral solver to try to solve this conformation
# solver will rotate each dihedral to reach a conformation where no atoms on one side of the dihedral are within a specified distance of any atom on the other side of the dihedral

# by default, this only looks at atoms from residues with resid up to the current monomer
# we are using a cutoff of 0.7 angrstroms
# this means that when you attempt to solve dihedral in monomer 5, it splits monomer 5 into two halves, separated by the dihedral
# and it tries to find a conformation where two things are true:
#  no atom on one side of the dihedral is within 0.7 A of any atom on the other side of the dihedral
#  no atom in monomer 5 is within 0.7 A of any atom on monomers with resid between zero and 5

# it does not look beyond residue 5.  
# the reason for this is that we want to solve the polymer starting from one end, and working along the chain

cutoff = 0.7 # we the generated conformation to have no two atoms within 0.7 A of each other
# this sort of distance is typically far enough to remove very high energy from overlapping atoms so that the structure can be further refined with energy minimzation

polymer4.dihedral_solver(CN_dihedrals,dummies='CX,NX',cutoff=cutoff)  #TODO why have I got a comma in there
# this will iterate through every bond in C2-N1 bond, rotating them so that atoms one either side of the bond are always more than 0.7 A apart.
# when it checks interatom distances, it ignores dummy atoms
# we have to specify the names of dummy atoms

Saver = PDB(polymer4)

Saver.cleanup() # center in box
Saver.save(dummyAtoms='CX NX',fname='polymer_04_solved') # save, excluding dummy atoms

# When you examine polymer four, you can see that the resulting structure is much more reasonable than polymer 1.  
# The polymer is arranged in a twisting helix.  The loops are much less tightly packed, and the distances between adjacent loops are much more reasonable.
# Furthermore, because we did not use a linear extend, the bonds and angles in this structure are consistent with those in the original monomers
# However, dihedal_solver() is algorithmic, and so it has still resulted in a highly ordered structure.  
# This polymer may be reasonable for illustrative purposes, but we may wish to generate an ensemble of less ordered structure for production simulations.

###############################################################################
#                                                                             #
#  Polymer five: Generating a randomised starting conformation shuffler()     #
#                                                                             #
###############################################################################

# in this example, we will make a copy of polymer 1, then use shuffler() to generate a random conformation
# we will then use dihedral_solver() to convert the random conformation into something reasonable

#first, make a copy of polymer one

polymer5=polymer1.copy()


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

NC_dihedrals=polymer5.gen_pairlist(J='N1',K='C1',first_resid=1,last_resid=127,mult=6,same_res=False) 
# these dihedrals are centered on the bond that links two monomers, so we need to set same_res to False
# we have also set the multiplicity of these dihedrals to 6, to increase the search space and improve the chance of finding a valid conformation

alkane_dhedrals=polymer5.gen_pairlist(J='C1',K='C2',first_resid=1,last_resid=128,mult=3)

polymer5.shuffler(NC_dihedrals,)  # this generates a random conformation 

# lets look at the intermediate conformation we got by random shuffling

Saver = PDB(polymer5)
Saver.cleanup() # center in box
Saver.save(dummyAtoms='CX NX',fname='polymer_05_shuffled') # save, excluding dummy atoms

# now solvce the conformation
polymer5.dihedral_solver(alkane_dhedrals,dummies='CX NX',cutoff=1) # this converts the shuffled conformation into one without overlapping atoms

Saver = PDB(polymer5)
Saver.cleanup() # center in box
Saver.save(dummyAtoms='CX NX',fname='polymer_05_shuffled_then_solved') # save, excluding dummy atoms

# Polymer five has two output files, shuffled, and shuffled_then_solved. 
# When you examine polymer five, you can see that the two conformations are extended random coils, and is quite different to the conformations of polymers one to four.

# The shuffled conformation has some steric clashes.  For example, monomers 47 and 54 are overlapping.
# The shuffled_then_solved conformation has been adjusted to remove these clashes.

# TODO make branched tutorial