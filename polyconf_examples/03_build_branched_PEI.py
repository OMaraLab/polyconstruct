#!/usr/bin/env python

from polyconf.Monomer import Monomer
from polyconf.Polymer import Polymer
from polyconf.PDB import PDB

import random

random.seed(10) # by setting the random seed we can ensure the random conformations are replicable.  

###############################################################################
#                                                                             #
#  Branched PEI                                                               #
#                                                                             #
###############################################################################

# In this example, we will follow the process used to generate a linear PEI, 128 monomers long

# but then we will add branches every 5 monomers 

PEI_branched=Polymer(Monomer('PEI_start.pdb')) 

adds=126 

for i in range (0,adds):
    PEI_branched.extend( 
            Monomer('PEI_monomer.pdb'),
            n=PEI_branched.maxresid(), 
            nn=PEI_branched.newresid(),
            names=dict(P='C1',Q='CX',R='NX',S='N1'), 
            joins=[('N1','C1')],
            ) 

# add the final monomer

PEI_branched.extend(
            Monomer('PEI_end.pdb'),
            n=PEI_branched.maxresid(), 
            nn=PEI_branched.newresid(),
            names=dict(P='C1',Q='CX',R='NX',S='N1'),
            joins=[('N1','C1')],
            )


# now we have a linear chain, we can add our branches
# for this example each branch will be three units long

branch_bonds=[] # make a list to keep track of all of our branches
branch=1

for n in range(5,128,5):

    branch_resid=PEI_branched.newresid()
    branch_bonds+=[(n,branch_resid)]
    # add the first monomer of the branch
    PEI_branched.extend( 
            Monomer('PEI_monomer.pdb'),
            n=n, 
            nn=branch_resid,
            names=dict(P='C1',Q='H1',R='NX',S='N1'),
            joins=[('N1','C1')],
            beta=branch # we are using beta factors to label our branches
            ) 

    PEI_branched.renamer(n,'H1')

    # add the second monomer of the branch

    PEI_branched.extend(
            Monomer('PEI_monomer.pdb'),
            n=PEI_branched.maxresid(), 
            nn=PEI_branched.newresid(),
            names=dict(P='C1',Q='CX',R='NX',S='N1'),
            joins=[('N1','C1')],
            beta=branch

            )

    # add the third monomer of the branch

    PEI_branched.extend(
            Monomer('PEI_end.pdb'),
            n=PEI_branched.maxresid(), 
            nn=PEI_branched.newresid(),
            names=dict(P='C1',Q='CX',R='NX',S='N1'),
            joins=[('N1','C1')],
            beta=branch

            )

    branch+=1 # change the branch number for the next branch

# Generate lists of dihedrals 

mainchain_NC= PEI_branched.gen_pairlist(J='N1',K='C1',first_resid=1,last_resid=127,mult=3,same_res=False) 
mainchain_alkanes= PEI_branched.gen_pairlist(J='C1',K='C2',first_resid=1,last_resid=128,mult=3)

# we are going to manually generate a dihedral list for the connections between the main chain and the branches
branch_dihedrals=[{'J': 'N1', 'J_resid': i, 'K': 'C1', 'K_resid': j, 'mult': 3} for i,j in branch_bonds]

branch_alkanes=PEI_branched.gen_pairlist(J='C1',K='C2',first_resid=129,last_resid=PEI_branched.maxresid(),mult=3)

for conf in range(1,6): 

    print(f'generating conformation {conf}')
    newconf=PEI_branched.copy() # make a duplicate of the original polymer

    print(f'shuffling backbone')

    newconf.shuffler(mainchain_NC)  # as before, randomise the NC dihedrals to generate a random conformation 

    print(f'solving backbone')

    newconf.dihedral_solver(mainchain_alkanes,dummies='CX NX X*',cutoff=0.9) # solve by rotating the alkane dihedrals

    print(f'shuffling branchpoints')

    newconf.shuffler(branch_dihedrals) 

    print(f'solving branchpoints')

    newconf.dihedral_solver(branch_dihedrals,dummies='CX NX X*',cutoff=0.9,backwards_only=False) 

    print(f'solving branches')

    newconf.dihedral_solver(branch_alkanes,dummies='CX NX X*',cutoff=0.9,backwards_only=False) 



    Saver = PDB(newconf)
    Saver.cleanup()
    Saver.save(dummies='CX NX X1',fname=f'PEI_branched_conformation_{conf}') 

# This should produce an ensemble of five starting conformations

# when you visualize these, try colouring by beta factor
# you can see each branch chain is labelled with different beta factor 
# this can be useful for manipulating and keeping track of branches