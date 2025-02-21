#!/usr/bin/env python

# Construction of an ethylamine dendrimer

# Import required Classes from PolyTop
from polytop.Junction import Junction
from polytop.Monomer import Monomer
from polytop.Visualize import Visualize
from polytop.Polymer import Polymer
from polytop.Topology import Topology


# ----- Load in monomer Topologies from the extended monomer ITP files -----

# Topology format is 'gromos' by default, but it is recommended to specify the
# format for clarity and readability.
core_mono = Topology.from_ITP("data/dendrimer_core.itp", format="gromos")
bifurcating_mono = Topology.from_ITP("data/dendrimer_bifurcating.itp", format="gromos")
terminal_mono = Topology.from_ITP("data/dendrimer_terminal.itp", format="gromos")


# ----- Create Junctions for the different Monomers to join to and from -----

# Provide the bonding atom and the leaving atom, in that order, for the Junction
# The two atoms used to create a Junction MUST have a bond between them.

# Note that junctions are specified for extend() by their name attribute and
# NOT by the variable name they are assigned to (which are instead used to pass
# them into a Monomer object).

# You will see that all of the Junctions created below have a unique name to
# prevent randomisation or uncertainty in the extension, and thus ensure repeatability.

# Create junctions for different 'levels' of monomers, depending on what level
# of branching they exist at in the dendrimer (e.g. central, first branching
# level, second branching level or terminal)

# Junctions for the 'central' Monomer that will initiate the Polymer
# One junction for each of the branches radiating out
central1 = Junction(core_mono.get_atom("N2"), core_mono.get_atom("C9"), name="C1")
central2 = Junction(core_mono.get_atom("N2"), core_mono.get_atom("C7"), name="C2")
central3 = Junction(core_mono.get_atom("N1"), core_mono.get_atom("C2"), name="C3")
central4 = Junction(core_mono.get_atom("N1"), core_mono.get_atom("C3"), name="C4")

# Junctions for the first level of branching with a 'bifurcating' Monomer
# The 'to' Junction will join to one of the 'central' Monomer Junctions, while
# the 'from1' and 'from2' Junctions will join to incoming second level of branching Monomers
b1 = Junction(bifurcating_mono.get_atom("C6"), bifurcating_mono.get_atom("N2"), name = "to")
b2a = Junction(bifurcating_mono.get_atom("N1"), bifurcating_mono.get_atom("C2"), name = "from1")
b2b = Junction(bifurcating_mono.get_atom("N1"), bifurcating_mono.get_atom("C3"), name = "from2")

# Junctions for the second level of branching with a 'bifurcating' Monomer
# The 'to2' Junction will join to one of the first level branching 'bifurcating'
# Monomer Junctions (i.e. from1 or from2), while the 'from12' and 'from22'
# Junctions will join to incoming terminal Monomers
c1 = Junction(bifurcating_mono.get_atom("C6"), bifurcating_mono.get_atom("N2"), name = "to2")
c2a = Junction(bifurcating_mono.get_atom("N1"), bifurcating_mono.get_atom("C2"), name = "from12")
c2b = Junction(bifurcating_mono.get_atom("N1"), bifurcating_mono.get_atom("C3"), name = "from22")

# Junctions for the 'terminal' Monomer that will terminate the Polymer
# Note how only 1 Junction is created for the terminal Monomer, as this will
# act to cap the Polymer and prevent further extension as there will be no
# remaining Junctions.
# The 'term' Junction will join to one of the second level branching 'bifurcating'
# Monomer Junctions (i.e. from12 or from22)
t = Junction(terminal_mono.get_atom("C1"), terminal_mono.get_atom("N1"), name = "term")


# ----- Create Monomers from their Topologies and any specified Junctions -----

# Note that different 'levels' of monomers have different names for their junctions.
# This ensures that layers are added on sequentially. 
# If all monomers have the same junction names the polymerisation will be random
# and often defaults to a linear shape!
central = Monomer(core_mono, [central1, central2, central3, central4])
bifur1 = Monomer(bifurcating_mono, [b1, b2a, b2b])
bifur2 = Monomer(bifurcating_mono, [c1, c2a, c2b])
cap = Monomer(terminal_mono, [t])


# ----- Start the Polymer with one Monomer -----

# Only a Polymer object has access to the 'extend()' and 'extra_bond()'
# functions needed to join monomers together, so the first Monomer object must
# be converted to a Polymer object to build from.
# Since the polymer will be built starting from the 'central' Monomer and
# radiating outwards, the 'central' Monomer will be used to start the Polymer.
polymer = Polymer(central)


# ----- Extending the Polymer to reach the desired structure -----

# The extend() function joins an incoming Monomer object to the Polymer object
# that is invoking the method (in this case 'polymer').
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

# Extend first layer of dendrimer
polymer.extend(bifur1, from_junction_name="C1", to_junction_name="to")
polymer.extend(bifur1, from_junction_name="C2", to_junction_name="to")
polymer.extend(bifur1, from_junction_name="C3", to_junction_name="to")
polymer.extend(bifur1, from_junction_name="C4", to_junction_name="to")

# Extend second layer of dendrimer
# Note how the first and second levels of branching have different Monomers
# with different Junction names to ensure the correct structure is always
# formed. Junction and Monomer name ambiguity or redundancy *WILL* cause
# formation of random, unreplicatable polymer topologies.
# These extends could be completed within a for loop to reduce the amount of
# code, but have been expanded for this example to make sure that all commands
# are clear and explicit.
polymer.extend(bifur2, from_junction_name="from1", to_junction_name="to2")
polymer.extend(bifur2, from_junction_name="from2", to_junction_name="to2")
polymer.extend(bifur2, from_junction_name="from1", to_junction_name="to2")
polymer.extend(bifur2, from_junction_name="from2", to_junction_name="to2")
polymer.extend(bifur2, from_junction_name="from1", to_junction_name="to2")
polymer.extend(bifur2, from_junction_name="from2", to_junction_name="to2")
polymer.extend(bifur2, from_junction_name="from1", to_junction_name="to2")
polymer.extend(bifur2, from_junction_name="from2", to_junction_name="to2")

# Finish polymer by extending on the capping monomer
# Note how each level of extension radiating outward from the central Monomer
# has double the number of extend() calls as the branching doubles each layer.
polymer.extend(cap, from_junction_name="from12", to_junction_name="term")
polymer.extend(cap, from_junction_name="from22", to_junction_name="term")
polymer.extend(cap, from_junction_name="from12", to_junction_name="term")
polymer.extend(cap, from_junction_name="from22", to_junction_name="term")
polymer.extend(cap, from_junction_name="from12", to_junction_name="term")
polymer.extend(cap, from_junction_name="from22", to_junction_name="term")
polymer.extend(cap, from_junction_name="from12", to_junction_name="term")
polymer.extend(cap, from_junction_name="from22", to_junction_name="term")
polymer.extend(cap, from_junction_name="from12", to_junction_name="term")
polymer.extend(cap, from_junction_name="from22", to_junction_name="term")
polymer.extend(cap, from_junction_name="from12", to_junction_name="term")
polymer.extend(cap, from_junction_name="from22", to_junction_name="term")
polymer.extend(cap, from_junction_name="from12", to_junction_name="term")
polymer.extend(cap, from_junction_name="from22", to_junction_name="term")
polymer.extend(cap, from_junction_name="from12", to_junction_name="term")
polymer.extend(cap, from_junction_name="from22", to_junction_name="term")

# Optionally (but recommended), check the netcharge of the monomers is
# preserved (in this case, close to 0)
print(polymer.topology.netcharge)


# ----- Save the polymer dendrimer to an itp file -----

# Optionally, use Visualize to generate an image of the structure with RDKit
# for an easy visual structure check
polymer.save_to_file('data/ethylamine_dendrimer.json') # optional, text dump in a dictionary format
polymer.topology.to_ITP('data/ethylamine_dendrimer.itp') # write the itp to be used for simulation
Visualize.polymer(polymer,infer_bond_order=False).draw2D('data/ethylamine_dendrimer.png',(400,300)) # optional, visualize the structure in 2D