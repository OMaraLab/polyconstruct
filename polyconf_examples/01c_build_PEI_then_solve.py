#!/usr/bin/env python

from polyconf.Monomer import Monomer
from polyconf.Polymer import Polymer
from polyconf.PDB import PDB


###############################################################################
#                                                                             #
#  Polymer four: Building a reasonable conformation with dihedral_solver()    #
#                                                                             #
###############################################################################



# in this example, we will make a copy through the same process as in polymer 1, then use dihedral_solver() to generate a more reasonable starting conformation

polymer4=Polymer(Monomer('PEI_start.pdb')) 

adds=126 

for i in range (0,adds):
    polymer4.extend( 
            Monomer('PEI_monomer.pdb'), 
            n=polymer4.maxresid(),
            nn=polymer4.newresid(), 
            names=dict(P='C1',Q='CX',R='NX',S='N1'),
            joins=[('N1','C1')],
            ) 

polymer4.extend(
            Monomer('PEI_end.pdb'),
            n=polymer4.maxresid(), 
            nn=polymer4.newresid(),
            names=dict(P='C1',Q='CX',R='NX',S='N1'),
            joins=[('N1','C1')],
            )

# we need to choose a set of dihedrals to shuffle

# I am going to choose every dihedral centered on the torsion between atoms N1 and C1.  This is the connection between monomers.

# we use gen pairlist to generate a list of all of these dihedrals
# note that we only specify the two atoms in the center of the dihedral

dihedrals=polymer4.gen_pairlist(J='N1',K='C1',first_resid=1,last_resid=127,mult=3,same_res=False) 

# next, we use dihedral solver to try to solve this conformation
# solver will rotate each dihedral to attempt to find a conformation where no atoms on one side of the dihedral are within a specified distance of any atom on the other side of the dihedral

cutoff = 0.9

# we want the generated conformation to have no two atoms within 0.9 A of each other

# this means that when you attempt to solve dihedral between monomers 4 and 5, it splits the polymer into two halves, separated by the dihedral
# then rotates the torsion tries to find a conformation where no atom on one side of the dihedral is within 0.9 A of any atom on the other side of the dihedral, considering only atoms in monomers 1 to 5

# it does not look beyond residue 5.  
# the reason for this is that we want to solve the polymer starting from one end, and working along the chain

# In my experience, this sort of distance is typically sufficient to produce a structure suitable for energy minimization, followed by equilibration.

polymer4.dihedral_solver(dihedrals,dummies='CX NX',cutoff=cutoff) 
# this will iterate through every N1-C1 torsion, rotating them so that atoms one either side of the torsion are always more than 0.7 A apart.
# when it checks interatom distances, it ignores dummy atoms provided we have specified the names of dummy atoms

Saver = PDB(polymer4)

Saver.cleanup() # center in box
Saver.save(dummies='CX NX',fname='PEI_linear_04_solved') # save, excluding dummy atoms

# When you examine polymer four, you can see that the resulting structure is much more reasonable than polymer 1.  
# The polymer is arranged in a twisting helix.  The loops are much less tightly packed, and the distances between adjacent loops are much more reasonable.
# Furthermore, because we did not use a linear extend, the bonds and angles in this structure are consistent with those in the original monomers
# However, dihedal_solver() is algorithmic, and so it has still resulted in a highly ordered structure.  
# This polymer may be reasonable for illustrative purposes, but we may wish to generate an ensemble of less ordered structure for production simulations.
