#!/usr/bin/env python
 
# Construction of a simple linear homopolymer of PEI

# This polymer will be 20 monomers long and constructed by joining Atom N71 to
# C5, and discarding the atoms C51 and N7 (and any other atoms, such as
# hydrogens, connected to them but not the rest of the polymer).

# Import required Classes from PolyTop
from polytop.Junction import Junction
from polytop.Monomer import Monomer
from polytop.Visualize import Visualize
from polytop.Polymer import Polymer
from polytop.Topology import Topology


# ----- Load in the monomer Topology from the extended PEI monomer ITP file -----

# Topology format is 'gromos' by default, but it is recommended to specify the
# format for clarity and readability.
top = Topology.from_ITP("data/pei.itp", format="gromos")


# ----- Create a Junction to join 'to' and another to join 'from' -----

# Provide the bonding atom and the leaving atom, in that order, for the Junction
# The two atoms used to create a Junction MUST have a bond between them.
# Note that junctions are specified for extend() by their name attribute (in
# this example 'to' and 'from', and NOT by the variable name they are assigned
# to 'to_j' and 'from_j', which are used to pass them into a Monomer object)
to_j = Junction(top.get_atom("N71"), top.get_atom("C51"), name = "from")
from_j = Junction(top.get_atom("C5"), top.get_atom("N7"), name = "to")

# ----- Create a Monomer from the Topology and a list of the Junctions -----

# Note that you could also create an 'initial' Monomer with the 'from_j'
# Junction only and a 'terminal' Monomer with the 'to_j' only, but this is not
# necessary, as any Junctions which are not extended to or from with extend()
# will not lose their leaving atoms. However, this is a beneficial approach for
# more complex polymers, such as with branching, to prevent ambiguity during
# extension. For examples of this, see tutorial scripts in
# `dendrimer_ethlyamine.py` or `star_PEG.py`
monomer = Monomer(top, [to_j, from_j])

# ----- Start the Polymer with one Monomer -----

# Only a Polymer object has access to the 'extend()' and 'extra_bond()'
# functions needed to join monomers together, so the first Monomer object must
# be converted to a Polymer object to build from.
polymer = Polymer(monomer)


# ----- Extend the Polymer to the desired length (in this case 20) -----

# The extend() function joins an incoming Monomer object (in this case 'monomer')
# to the Polymer object that is invoking the method (in this case 'polymer').
# The argument 'from_junction_name' specifies which Junction (by its name) in
# the Polymer will form a new bond to the incoming Monomer, while the argument
# 'to_junction_name' is the Junction in the incoming Monomer that will form a
# bond to the Polymer. Such that the Polymer Junction's 'monomer atom' will
# join to the Monomer Junction's 'monomer atom' and both of their 'leaving atoms'
# will be discarded.
# To ensure that you generate a valid topology, the two Junctions you join
# should have corresponding atoms and the same bond, angle and dihedral
# parameters, which will be ensured by using an appropriately extended topology
# and defining Junctions correctly.
for i in range(19):
    polymer.extend(monomer, from_junction_name="from", to_junction_name="to")

# ----- Save the polymer to an itp file -----
# Optionally, use Visualize to generate an image of the structure with RDKit
# for an easy visual structure check 
polymer.topology.title = "pei polymer" # optional but good for identifying files, renames the ITP header and image
polymer.save_to_file('data/pei_linear_polymer.json') # optional, text dump in a dictionary format
polymer.topology.to_ITP('data/pei_linear_polymer.itp') # write the itp to be used for simulation
Visualize.polymer(polymer,infer_bond_order=False).draw2D('data/pei_linear_polymer.png',(400,300)) # optional, visualize the structure in 2D