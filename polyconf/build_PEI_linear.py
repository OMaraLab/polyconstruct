#!/usr/bin/env python

from polyconf.monomer import Monomer
from polyconf.polymer import Polymer
from polyconf.PDB import PDB

# polymer=Polymer(Monomer('PEI_start.pdb')) # initialise
# imax=127 # lay 127 additional monomers

# # first, try building with normal extend

# for i in range (1,imax+1):
#     if not i==imax: # do this for everything except the first one
#         polymer.extend( # extend with one monomer, aligned along this step's linearization vector
#             Monomer('PEI_monomer.pdb'), # extend with this monomer
#             i, # extend existing residue i
#             i+1, # incoming monomer will have resid i+1
#             names=dict(P1='CX',P2='C1',Q1='N1',Q2='NX'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
#             joins=[('N1','C1')],# new connection between N1_i and C1_i+1 
#             ) 
#     else:
#         # add final monomer
#         polymer.extend(Monomer('PEI_end.pdb'),i,i+1,
#             names=dict(P1='CX',P2='C1',Q1='N1',Q2='NX'),
#             joins=[('N1','C1')],)
# #    polymer.genconf([('C2','N1')],length=10, cutoff = 0.5)
#     'C1' 'NX' 'CX' 'N1'

# Saver = PDB(polymer)
# Saver.cleanup() # center in box
# Saver.save(dummyAtoms='CX NX',fname='polymer_01_vanilla-extend') # save, excluding dummy atoms


# # first, try building with linear extend


# polymer=Polymer(Monomer('PEI_start.pdb')) # initialise
# imax=127 # lay 127 additional monomers
# ortho=[1,0,0] # linearization vector



# for i in range (1,imax+1):
#     if not i==imax: # do this for everything except the first one
#         polymer.extend(
#             Monomer('PEI_monomer.pdb'), # extend with this monomer
#             i, # extend existing residue i
#             i+1, # incoming monomer will have resid i+1
#             names=dict(P1='CX',P2='C1',Q1='N1',Q2='NX', # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
#             R='C1',S='N1'), # monomer vector is vector C1_i+1 N1_i+1,
#             linearise=True, # linearize incoming monomer so monomer vector lies along along ortho
#             ortho=ortho,# specify linearization vector
#             joins=[('N1','C1')],# new connection between N1_i and C1_i+1 
#             ) 
#     else:
#         # add final monomer
#         polymer.extend(Monomer('PEI_end.pdb'),i,i+1,
#             names=dict(P1='CX',P2='C1',Q1='N1',Q2='NX'),
#             joins=[('N1','C1')],)
# #    polymer.genconf([('C2','N1')],length=10, cutoff = 0.5)
#     'C1' 'NX' 'CX' 'N1'

# Saver = PDB(polymer)
# Saver.cleanup() # center in box
# Saver.save(dummyAtoms='CX NX',fname='polymer_02_simple-linear-extend') # save, excluding dummy atoms


# third, try building with linear extend, but alternating between two ortho vectors


polymer=Polymer(Monomer('PEI_start.pdb')) # initialise
imax=23 # lay 24 additional monomers
alternator=True # lay monomers along alternating vectors to get a clean linear conformation, because shuffle isn't working



for i in range (1,imax+1):
    if not i==imax: # do this for everything except the first one
        if alternator:
            ortho=[1,0,0] # linearization vector A
        else:
            ortho=[0,1,0] # linearization vector B
        alternator=not alternator # flip the alternator

        polymer.extend(Monomer( # extend with one monomer, aligned along this step's linearization vector
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
        polymer.extend(Monomer('PEI_end.pdb'),i,i+1,
            names=dict(P1='CX',P2='C1',Q1='N1',Q2='NX'),
            joins=[('N1','C1')],)
#    polymer.genconf([('C2','N1')],length=10, cutoff = 0.5)
    'C1' 'NX' 'CX' 'N1'

Saver = PDB(polymer)
Saver.cleanup() # center in box
Saver.save(dummyAtoms='CX NX',fname='polymer_03_alternating-linear-extend') # save, excluding dummy atoms


alkanes=polymer.gen_pairlist(a1='C1',a2='C2',first_resid=1,last_resid=24,mult=6)

print(alkanes[0:5])
polymer.shuffler(alkanes,dummy='CX NX',cutoff=0.5)
Saver = PDB(polymer)
Saver.cleanup() # center in box
Saver.save(dummyAtoms='CX NX',fname='polymer_04_shuffled') # save, excluding dummy atoms


polymer.dihedral_solver(alkanes,dummy='CX,NX',cutoff=0.5)

Saver = PDB(polymer)
Saver.cleanup() # center in box
Saver.save(dummyAtoms='CX NX',fname='polymer_05_solved') # save, excluding dummy atoms
