#!/usr/bin/env python

from polyconf.Monomer import Monomer
from polyconf.Polymer import Polymer
from polyconf.PDB import PDB
import random
random.seed(0) 

print ('creating an atactic PMMA polymer using a linear extend()') 

dummies="CMA CN CP CQ"

Monomer_D='MMAD_bonds.pdb'
Monomer_L='MMAL_bonds.pdb'

composition=25*[Monomer_D] + 25*[Monomer_L] # a list containing 50 monomers, 25 with each taticity
composition = random.sample(composition,len(composition)) # randomise the order of the monomers

PMMA_atactic_linear=Polymer(Monomer(composition[0])) # initialise with first monomer

for monomer in composition[1:]: # extend with the remain 49 monomers
        PMMA_atactic_linear.extend( # extend with one monomer, aligned along this step's linearization vector
            Monomer(monomer), # extend with this monomer
            n=PMMA_atactic_linear.maxresid(), # extend existing residue i
            nn=PMMA_atactic_linear.newresid(), # incoming monomer will have resid i+1
            names=dict(Q='CA',P='CMA',S='C',R='CN',V1='CN',V2='C'), # C1_i+1 fit to CX_i,
            joins=[('C','CA')],# new connection between N1_i and C1_i+1 
            linearise=True,
            ortho=[1,0,0] 
            ) 

Saver = PDB(PMMA_atactic_linear)
Saver.cleanup() 
Saver.save(dummies=dummies,fname='PMMA_atactic_linear-extend') # save, excluding dummies atoms

# when you look at this polymer, you can see the alkane backbone is extremely regular
# while this is not realistic, we can use it to generate reasonable conformations 

alkanes=PMMA_atactic_linear.gen_pairlist(J='CA',K='C',first_resid=1,same_res=True,last_resid=50,mult=6)
sidechains=PMMA_atactic_linear.gen_pairlist(J='CA',K='CB',first_resid=1,same_res=True,last_resid=50,mult=6)

dh=[] #combine alkanes and sidechains into one list of dihedrals, so we can shuffle the alkane and dihedral in each monomer in order

for i in range(0,50):
    dh += [alkanes[i]]
    dh += [sidechains[i]]

print ('solving pseudolinear atactic PMMA ') 

PMMA_atactic_linear.dihedral_solver(dh,dummies=dummies,cutoff=1.1)

# unsurprisingly it was very easy to solve this pseudolinear conformation

Saver = PDB(PMMA_atactic_linear)
Saver.cleanup() 
Saver.save(dummies=dummies,fname='PMMA_atactic_linear-solved') # save, excluding dummies atoms

print ('shuffle attempt 01') 

PMMA_atactic_linear_01 = PMMA_atactic_linear.copy()

PMMA_atactic_linear_01.shuffler(dh)
PMMA_atactic_linear_01.dihedral_solver(dh,dummies=dummies,cutoff=1.1)

Saver = PDB(PMMA_atactic_linear_01)
Saver.cleanup() 
Saver.save(dummies=dummies,fname='PMMA_atactic_linear-shuffled_and_solved_01') # save, excluding dummies atoms

# this generated a valid conformation
# lets repeat the process three more times

print ('shuffle attempt 02') 

PMMA_atactic_linear_02 = PMMA_atactic_linear.copy()

PMMA_atactic_linear_02.shuffler(dh)
PMMA_atactic_linear_02.dihedral_solver(dh,dummies=dummies,cutoff=1.1)

Saver = PDB(PMMA_atactic_linear_02)
Saver.cleanup() 
Saver.save(dummies=dummies,fname='PMMA_atactic_linear-shuffled_and_solved_02') # save, excluding dummies atoms


print ('shuffle attempt 03') 

PMMA_atactic_linear_03 = PMMA_atactic_linear.copy()

PMMA_atactic_linear_03.shuffler(dh)
PMMA_atactic_linear_03.dihedral_solver(dh,dummies=dummies,cutoff=1.1)

Saver = PDB(PMMA_atactic_linear_03)
Saver.cleanup() 
Saver.save(dummies=dummies,fname='PMMA_atactic_linear-shuffled_and_solved_03') # save, excluding dummies atoms

print ('shuffle attempt 04') 

PMMA_atactic_linear_04 = PMMA_atactic_linear.copy()

PMMA_atactic_linear_04.shuffler(dh)
PMMA_atactic_linear_04.dihedral_solver(dh,dummies=dummies,cutoff=1.1)

Saver = PDB(PMMA_atactic_linear_04)
Saver.cleanup() 
Saver.save(dummies=dummies,fname='PMMA_atactic_linear-shuffled_and_solved_04') # save, excluding dummies atoms

# all four attempts were able to be solved, generating four possible starting conformations for energy minimzation
