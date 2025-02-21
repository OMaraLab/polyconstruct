#!/usr/bin/env python

from polyconf.Monomer import Monomer
from polyconf.Polymer import Polymer
from polyconf.PDB import PDB
import random
random.seed(10) 

###############################################################################
#                                                                             #
#  Polymer three: Building an atactic PMMA polymer                            #
#                                                                             #
###############################################################################

print ('creating an atactic PMMA polymer') 

dummies="CMA CN CP CQ"
Monomer_D='MMAD_bonds.pdb'
Monomer_L='MMAL_bonds.pdb'

composition=25*[Monomer_D] + 25*[Monomer_L] # a list containing 50 monomers, 25 with each taticity
composition = random.sample(composition,len(composition)) # randomise the order of the monomers

PMMA_atactic=Polymer(Monomer(composition[0])) # initialise with first monomer

for monomer in composition[1:]: # extend with the remain 49 monomers
        PMMA_atactic.extend( # extend with one monomer, aligned along this step's linearization vector
            Monomer(monomer), # extend with this monomer
            n=PMMA_atactic.maxresid(), # extend existing residue i
            nn=PMMA_atactic.newresid(), # incoming monomer will have resid i+1
            names=dict(Q='CA',P='CMA',S='C',R='CN'), # C1_i+1 fit to CX_i,
            joins=[('C','CA')],# new connection between N1_i and C1_i+1 
            ) 

Saver = PDB(PMMA_atactic)
Saver.cleanup() # center in box
Saver.save(dummies=dummies,fname='PMMA_atactic_vanilla-extend') # save, excluding dummies atoms

# have a look at this initial polymer, compared to the vanilla isotactic and syndiotactic polymers
# because of the particular conformations of the two monomers, and the random order of the monomer composition
# the backbone of this polymer is much more coiled.absyou might think this is desifreable, as it reflects a less ordered starting conformation
# however, many monomers have overlapping atoms, and it is unlikely that we could minimise this starting confromation

alkanes=PMMA_atactic.gen_pairlist(J='CA',K='C',first_resid=1,same_res=True,last_resid=50,mult=6)
sidechains=PMMA_atactic.gen_pairlist(J='CA',K='CB',first_resid=1,same_res=True,last_resid=50,mult=6)

dh=[] #combine alkanes and sidechains into one list of dihedrals, so we can shuffle the alkane and dihedral in each monomer in order

for i in range(0,50):
    dh += [alkanes[i]]
    dh += [sidechains[i]]

print('attempting to solve initial conformation')
print('the random seed in this example script has been chosen specifically so that this process fails')
print('this is expected and is part of the tutorial')

PMMA_atactic.dihedral_solver(dh,dummies=dummies,cutoff=1.1) 

# this failed due to overlapping atoms

# we can attempt to rotate and shuffle dihedrals to find a valid conformation
# however the brute force search of dihedral space used by dihedral_solver() will fail for this particular random
# it is likely that repeated applications of shuffle and dihedral_solver will eventually find a reasonable conformation
# but we cannot know how many attempts this will take

print('shuffle attempt 01')

PMMA_atactic.shuffler(dh) 

PMMA_atactic.dihedral_solver(dh,dummies=dummies,cutoff=1.1 )

Saver = PDB(PMMA_atactic)
Saver.cleanup() # center in box
Saver.save(dummies=dummies,fname='PMMA_atactic_shuffle_01') # save, excluding dummies atoms

print('shuffle attempt 02')

PMMA_atactic.shuffler(dh) 
PMMA_atactic.dihedral_solver(dh,dummies=dummies,cutoff=1.1) 

Saver = PDB(PMMA_atactic)
Saver.cleanup() # center in box
Saver.save(dummies=dummies,fname='PMMA_atactic_shuffle_02') # save, excluding dummies atoms

print('shuffle attempt 03')

PMMA_atactic.shuffler(dh) 

PMMA_atactic.dihedral_solver(dh,dummies=dummies,cutoff=1.1) 

Saver = PDB(PMMA_atactic)
Saver.cleanup() # center in box
Saver.save(dummies=dummies,fname='PMMA_atactic_shuffle_03') # save, excluding dummies atoms


print('shuffle attempt 04')

PMMA_atactic.shuffler(dh) 

PMMA_atactic.dihedral_solver(dh,dummies=dummies,cutoff=1.1) 

Saver = PDB(PMMA_atactic)
Saver.cleanup() # center in box
Saver.save(dummies=dummies,fname='PMMA_atactic_shuffle_04') # save, excluding dummies atoms

# even after four steps, this is still unable to solve the conformation

# If your initial conformation contains many overlapping monomers, I have found you are much less likely to find valid starting conformations
# There are several steps we could take to resolve this
# One approach would be to reduce the cutoff, which allows conformations with smaller interatomic distances

# In this example, using a cutoff of 0.8 A will allow dihedral_solver() to generate a starting starting conformation after every step

# Reducing the cutoff means some atoms are likely to be unrealistically close, but an energy minimzation may address this and result in structures that can be used in molecular dynamics

# Another approach is to vary the multiplicity of your dihedrals, to alter the search parameters through dihedral space
# In this case the alkane dihedrals have been permitted to explore a multiplicity of 6, but changing it to 3 (restricting the search space) or to 12 (increasing the search space) may facilitate a solution

# You could also try solving the alkanes first, and then the sidechains, or doing the reverse

# By varying your approach, you will eventually stumble upon a method that identifies valid starting conformations

# However, I have found it is quite reliable to start by generating a psuedolinear conformation, and then shuffling and solving that conformation
# An example of this process is given in build_PMMA_04-atactic-linearized