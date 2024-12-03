#!/usr/bin/env python

from polyconf.monomer import Monomer
from polyconf.polymer import Polymer
from polyconf.PDB import PDB
import random

###############################################################################
#                                                                             #
#  Polymer one: Building an isotactic PMMA polymer                            #
#                                                                             #
###############################################################################

dummies="CMA CN CP CQ"

PMMA_isotactic=Polymer(Monomer('MMAD_bonds.pdb')) # initialise

imax=49 # we will lay 49 additional monomers

for i in range (0,imax):
        PMMA_isotactic.extend( # extend with one monomer, aligned along this step's linearization vector
            Monomer('MMAD_bonds.pdb'), # extend with this monomer
            n=PMMA_isotactic.maxresid(), # extend existing residue i
            nn=PMMA_isotactic.newresid(), # incoming monomer will have resid i+1
            names=dict(P1='CA',P2='CMA',Q1='C',Q2='CN',R='CN',S='C'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
            joins=[('C','CA')],# new connection between N1_i and C1_i+1 
            #linearise=True,
            #ortho=[1,0,0]
            ) 

alkanes=PMMA_isotactic.gen_pairlist(a1='C',a2='CA',first_resid=1,same_res=False,last_resid=49,mult=3)

Saver = PDB(PMMA_isotactic)
Saver.cleanup() # center in box
Saver.save(dummyAtoms=dummies,fname='PMMA_isotactic_vanilla-extend') # save, excluding dummy atoms

PMMA_isotactic.dihedral_solver(alkanes,dummy=dummies,cutoff=0.8) # this converts the shuffled conformation into one without overlapping atoms

Saver = PDB(PMMA_isotactic)
Saver.cleanup() # center in box
Saver.save(dummyAtoms=dummies,fname='PMMA_isotactic_vanilla-solved') # save, excluding dummy atoms

PMMA_isotactic.shuffler(alkanes)
PMMA_isotactic.dihedral_solver(alkanes,dummy=dummies,cutoff=1.1) # this converts the shuffled conformation into one without overlapping atoms

Saver = PDB(PMMA_isotactic)
Saver.cleanup() # center in box
Saver.save(dummyAtoms=dummies,fname='PMMA_isotactic_shuffled_and_solved') # save, excluding dummy atoms



###############################################################################
#                                                                             #
#  Polymer two: Building a syndiotactic PMMA polymer with                     #
#                                                                             #
###############################################################################

dummies="CMA CN CP CQ"

Monomer_D='MMAD_bonds.pdb'
Monomer_L='MMAL_bonds.pdb'
alternator=True

PMMA_syndiotactic=Polymer(Monomer(Monomer_D)) # initialise

imax=49 # we will lay 49 additional monomers

for i in range (0,imax):
        if alternator:
            monomer=Monomer_L
        else:
            monomer=Monomer_D
        alternator=not alternator
        PMMA_syndiotactic.extend( # extend with one monomer, aligned along this step's linearization vector
            Monomer(monomer), # extend with this monomer
            n=PMMA_syndiotactic.maxresid(), # extend existing residue i
            nn=PMMA_syndiotactic.newresid(), # incoming monomer will have resid i+1
            names=dict(P1='CA',P2='CMA',Q1='C',Q2='CN',R='CN',S='C'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
            joins=[('C','CA')],# new connection between N1_i and C1_i+1 
            #linearise=True,
            #ortho=[1,0,0]
            ) 

Saver = PDB(PMMA_syndiotactic)
Saver.cleanup() # center in box
Saver.save(dummyAtoms=dummies,fname='PMMA_syndiotactic_vanilla-extend') # save, excluding dummy atoms

alkanes=PMMA_syndiotactic.gen_pairlist(a1='C',a2='CA',first_resid=1,same_res=False,last_resid=49,mult=3)
PMMA_syndiotactic.dihedral_solver(alkanes,dummy=dummies,cutoff=1.1) # this converts the shuffled conformation into one without overlapping atoms

Saver = PDB(PMMA_syndiotactic)
Saver.cleanup() # center in box
Saver.save(dummyAtoms=dummies,fname='PMMA_syndiotactic_vanilla-solved') # save, excluding dummy atoms

sidechains=PMMA_syndiotactic.gen_pairlist(a1='CA',a2='CB',first_resid=1,same_res=True,last_resid=50,mult=6)

PMMA_syndiotactic.shuffler(alkanes)
PMMA_syndiotactic.dihedral_solver(alkanes,dummy=dummies,cutoff=1.1) # this converts the shuffled conformation into one without overlapping atoms
PMMA_syndiotactic.shuffler(sidechains)
PMMA_syndiotactic.dihedral_solver(sidechains,dummy=dummies,cutoff=1.1) # this converts the shuffled conformation into one without overlapping atoms

Saver = PDB(PMMA_syndiotactic)
Saver.cleanup() # center in box
Saver.save(dummyAtoms=dummies,fname='PMMA_syndiotactic_shuffled_and_solved') # save, excluding dummy atoms




###############################################################################
#                                                                             #
#  Polymer two: Building an atactic PMMA polymer                              #
#                                                                             #
###############################################################################

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
            names=dict(P1='CA',P2='CMA',Q1='C',Q2='CN',R='CN',S='C'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
            joins=[('C','CA')],# new connection between N1_i and C1_i+1 
            #linearise=True,
            #ortho=[0.5,-1,0]
            ) 

Saver = PDB(PMMA_atactic)
Saver.cleanup() # center in box
Saver.save(dummyAtoms=dummies,fname='PMMA_atactic_vanilla-extend') # save, excluding dummy atoms

alkanes=PMMA_atactic.gen_pairlist(a1='CA',a2='C',first_resid=1,same_res=True,last_resid=49,mult=12)
PMMA_atactic.dihedral_solver(alkanes,dummy=dummies,cutoff=1.1) # this converts the shuffled conformation into one without overlapping atoms

Saver = PDB(PMMA_atactic)
Saver.cleanup() # center in box
Saver.save(dummyAtoms=dummies,fname='PMMA_atactic_vanilla-solved') # save, excluding dummy atoms

sidechains=PMMA_atactic.gen_pairlist(a1='CA',a2='CB',first_resid=1,same_res=True,last_resid=50,mult=6)

PMMA_atactic.shuffler(alkanes)
PMMA_atactic.dihedral_solver(alkanes,dummy=dummies,cutoff=1.1) # this converts the shuffled conformation into one without overlapping atoms
# PMMA_atactic.shuffler(sidechains)
# PMMA_atactic.dihedral_solver(sidechains,dummy=dummies,cutoff=1.1) # this converts the shuffled conformation into one without overlapping atoms

Saver = PDB(PMMA_atactic)
Saver.cleanup() # center in box
Saver.save(dummyAtoms=dummies,fname='PMMA_atactic_shuffled_and_solved_01') # save, excluding dummy atoms

# PMMA_atactic=Polymer(Monomer(composition[0])) # initialise with first monomer

# for monomer in composition[1:]: # extend with the remain 49 monomers
#         PMMA_atactic.extend( # extend with one monomer, aligned along this step's linearization vector
#             Monomer(monomer), # extend with this monomer
#             n=PMMA_atactic.maxresid(), # extend existing residue i
#             nn=PMMA_atactic.newresid(), # incoming monomer will have resid i+1
#             names=dict(P1='CA',P2='CMA',Q1='C',Q2='CN',R='CN',S='C'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
#             joins=[('C','CA')],# new connection between N1_i and C1_i+1 
#             linearise=True,
#             ortho=[0.5,-1,0]
#             ) 

# PMMA_atactic.shuffler(alkanes)
# PMMA_atactic.dihedral_solver(alkanes,dummy=dummies,cutoff=1.1) # this converts the shuffled conformation into one without overlapping atoms

# Saver = PDB(PMMA_atactic)
# Saver.cleanup() # center in box
# Saver.save(dummyAtoms=dummies,fname='PMMA_atactic_shuffled_and_solved_02') # save, excluding dummy atoms

# PMMA_atactic=Polymer(Monomer(composition[0])) # initialise with first monomer

# for monomer in composition[1:]: # extend with the remain 49 monomers
#         PMMA_atactic.extend( # extend with one monomer, aligned along this step's linearization vector
#             Monomer(monomer), # extend with this monomer
#             n=PMMA_atactic.maxresid(), # extend existing residue i
#             nn=PMMA_atactic.newresid(), # incoming monomer will have resid i+1
#             names=dict(P1='CA',P2='CMA',Q1='C',Q2='CN',R='CN',S='C'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
#             joins=[('C','CA')],# new connection between N1_i and C1_i+1 
#             linearise=True,
#             ortho=[0.5,-1,0]
#             ) 
# PMMA_atactic.shuffler(alkanes)
# PMMA_atactic.dihedral_solver(alkanes,dummy=dummies,cutoff=1.1) # this converts the shuffled conformation into one without overlapping atoms

# Saver = PDB(PMMA_atactic)
# Saver.cleanup() # center in box
# Saver.save(dummyAtoms=dummies,fname='PMMA_atactic_shuffled_and_solved_03') # save, excluding dummy atoms
# PMMA_atactic=Polymer(Monomer(composition[0])) # initialise with first monomer

# for monomer in composition[1:]: # extend with the remain 49 monomers
#         PMMA_atactic.extend( # extend with one monomer, aligned along this step's linearization vector
#             Monomer(monomer), # extend with this monomer
#             n=PMMA_atactic.maxresid(), # extend existing residue i
#             nn=PMMA_atactic.newresid(), # incoming monomer will have resid i+1
#             names=dict(P1='CA',P2='CMA',Q1='C',Q2='CN',R='CN',S='C'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
#             joins=[('C','CA')],# new connection between N1_i and C1_i+1 
#             linearise=True,
#             ortho=[0.5,-1,0]
#             ) 
# PMMA_atactic.shuffler(alkanes)
# PMMA_atactic.dihedral_solver(alkanes,dummy=dummies,cutoff=1.1) # this converts the shuffled conformation into one without overlapping atoms
# PMMA_atactic=Polymer(Monomer(composition[0])) # initialise with first monomer

# for monomer in composition[1:]: # extend with the remain 49 monomers
#         PMMA_atactic.extend( # extend with one monomer, aligned along this step's linearization vector
#             Monomer(monomer), # extend with this monomer
#             n=PMMA_atactic.maxresid(), # extend existing residue i
#             nn=PMMA_atactic.newresid(), # incoming monomer will have resid i+1
#             names=dict(P1='CA',P2='CMA',Q1='C',Q2='CN',R='CN',S='C'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
#             joins=[('C','CA')],# new connection between N1_i and C1_i+1 
#             linearise=True,
#             ortho=[0.5,-1,0]
#             ) 
# Saver = PDB(PMMA_atactic)
# Saver.cleanup() # center in box
# Saver.save(dummyAtoms=dummies,fname='PMMA_atactic_shuffled_and_solved_04') # save, excluding dummy atoms

# PMMA_atactic.shuffler(alkanes)
# PMMA_atactic.dihedral_solver(alkanes,dummy=dummies,cutoff=1.1) # this converts the shuffled conformation into one without overlapping atoms

# Saver = PDB(PMMA_atactic)
# Saver.cleanup() # center in box
# Saver.save(dummyAtoms=dummies,fname='PMMA_atactic_shuffled_and_solved_05') # save, excluding dummy atoms
