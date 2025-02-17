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

# then, we will extend our polymer by adding new monomers
# we will repeat this 126 times

adds=126

for _ in range (0,adds): # repeat 126 times

    n=polymer1.maxresid()
    new_n=polymer1.newresid()

    polymer1.extend( # extend with one monomer
            Monomer('PEI_monomer.pdb'), # we are adding a new monomer using the coordinates from the file 'PEI_monomer.pdb'
            n=n,
            nn=new_n,
            names=dict(P='C1',Q='CX',R='NX',S='N1'), # this defines the mapping used to connect our monomers.  The new polymer will be created so the PR bond in the new monomer is colinear with the QS bond in the existing monomer
            joins=[('N1','C1')], # new connection between N1_n and C1_nn 
            ) 
    
    # after adding each monomer, we will flag dummy atoms
    # in this case, the dummy atoms are atoms CX in monomer n, and atom NX in monomer nn

    polymer1.renamer(n,'CX') # convert atom CX in monomer n to a dummy atom
    polymer1.renamer(new_n,'NX') # convert atom NX in monomer new_n to a dummy atom



# we now have a chain of 127 monomers
# now we add the final monomer
# the process is the same, but we are using a different pdb file

n = polymer1.maxresid()
new_n = polymer1.newresid()

polymer1.extend(
            Monomer('PEI_end.pdb'),
            n=polymer1.maxresid(), 
            nn=polymer1.newresid(),
            names=dict(P='C1',Q='CX',R='NX',S='N1'),
            joins=[('N1','C1')],
            )

polymer1.renamer(n,'CX')
polymer1.renamer(new_n,'NX')

# Now we can save our polymer as a pdb

# we start by making a PDB object, which will save our file
print('saving')

Saver = PDB(polymer1)

Saver.cleanup() # this puts our polymer in the center of the box
Saver.save(fname='PEI_linear_01_vanilla-extend') 

# When you examine polymer one , you can see that the resulting strucure is a tightly coiled helix.  
# This structure is highly ordered, and the turns of the helix are very close.
# It's not a very likely conformation

# In the next tutorial, we will use some other methods to make different conformations

