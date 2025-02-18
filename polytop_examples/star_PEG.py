# Construction of a 4-arm PEG star polymer from single monomeric units

# Import required Classes from PolyTop
from polytop.Junction import Junction
from polytop.Monomer import Monomer
from polytop.Visualize import Visualize
from polytop.Polymer import Polymer
from polytop.Topology import Topology


# ----- Load in monomer Topologies from the extended monomer ITP files -----

# Topology format is 'gromos' by default, but it is recommended to specify the
# format for clarity and readability.
ethanol = Topology.from_ITP("data/extended_ethanol.itp", format="gromos") # main arm monomer
methane = Topology.from_ITP("data/extended_methane.itp", format="gromos") # terminal monomer
neopentane = Topology.from_ITP("data/extended_neopentane.itp", format="gromos") # central monomer


# ----- Create Junctions for the different Monomers to join to and from -----

# Provide the bonding atom and the leaving atom, in that order, for the Junction
# The two atoms used to create a Junction MUST have a bond between them.

# Note that junctions are specified for extend() by their name attribute and
# NOT by the variable name they are assigned to (which are instead used to pass
# them into a Monomer object).

# You will see that all of the Junctions created below have a unique name to
# prevent randomisation or uncertainty in the extension, and thus ensure repeatability.

# Create unique Junctions for the Monomers that make up each of the 4 arms of
# the star, as well as the central Monomer and the terminal Monomers.

# Junctions for the central 'neopentane' Monomer that will initiate the Polymer
# One junction for each of the arms of the star
j1 = Junction(neopentane.get_atom("C1"), neopentane.get_atom("O1"), name = "branch1")
j2 = Junction(neopentane.get_atom("C3"), neopentane.get_atom("O2"), name = "branch2")
j3 = Junction(neopentane.get_atom("C4"), neopentane.get_atom("O3"), name = "branch3")
j4 = Junction(neopentane.get_atom("C5"), neopentane.get_atom("O4"), name = "branch4")

# Junctions for each of the ethanol 'arm' Monomers
# The 'oxy*' Junctions will join to one of the central 'neopentane' Monomer
# Junctions (branch*), while the 'carb*'Junctions will join to an incoming
# ethanol Monomer to extend that arm further. Such that branch1 Junction joins
# to oxy1 Junction, and the carb1 Junction of the same Monomer will join to
# another incoming Monomer's oxy1 Junction.
oxy_j1 = Junction(ethanol.get_atom("O1"), ethanol.get_atom("C1"), name = "oxy1")
carb_j1 = Junction(ethanol.get_atom("C3"), ethanol.get_atom("O2"), name = "carb1")
oxy_j2 = Junction(ethanol.get_atom("O1"), ethanol.get_atom("C1"), name = "oxy2")
carb_j2 = Junction(ethanol.get_atom("C3"), ethanol.get_atom("O2"), name = "carb2")
oxy_j3 = Junction(ethanol.get_atom("O1"), ethanol.get_atom("C1"), name = "oxy3")
carb_j3 = Junction(ethanol.get_atom("C3"), ethanol.get_atom("O2"), name = "carb3")
oxy_j4 = Junction(ethanol.get_atom("O1"), ethanol.get_atom("C1"), name = "oxy4")
carb_j4 = Junction(ethanol.get_atom("C3"), ethanol.get_atom("O2"), name = "carb4")

# Junctions for the 'terminal' Monomer that will terminate the Polymer's arms
# Note how only 1 Junction is created for the terminal Monomer, as this will
# act to cap the Polymer and prevent further extension as there will be no
# remaining Junctions.
# The 'term' Junction will join to the final 'ethanol' Monomer Junction (carb*)
term_j = Junction(methane.get_atom("C1"), methane.get_atom("O1"), name = "term")


# ----- Create Monomers from their Topologies and any specified Junctions -----

# Note that a unique Monomer has been created for each of the 4 arms with
# different Junctions, which ensures that only Monomers of the same type are
# joined to create the same arm, which means that all 4 arms will be the same
# desired length.
# If all Monomers had Junctions with the same names the polymerisation would be
# random and each of the arms would be different lengths. The lengths of each
# branch would also vary each time you built the Polymer, which is terrible for
# scientific replicability!
e1 = Monomer(ethanol, [oxy_j1, carb_j1])
e2 = Monomer(ethanol, [oxy_j2, carb_j2])
e3 = Monomer(ethanol, [oxy_j3, carb_j3])
e4 = Monomer(ethanol, [oxy_j4, carb_j4])

central = Monomer(neopentane, [j1, j2, j3, j4])

terminal = Monomer(methane, [term_j]) # only needs one junction to join to the ends of each arm


# ----- Start the Polymer with one Monomer -----

# Only a Polymer object has access to the 'extend()' and 'extra_bond()'
# functions needed to join monomers together, so the first Monomer object must
# be converted to a Polymer object to build from.
# Since the polymer will be built starting from the 'central' Monomer and
# radiating outwards, the 'central' Monomer will be used to start the Polymer.
four_polymer = Polymer(central)


# ----- Extending the Polymer to reach the desired structure -----

# The extend() function joins an incoming Monomer object to the Polymer object
# that is invoking the method (in this case 'four_polymer').
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

# attach three ethanols to each of the four junctions (branch1-branch4) of the
# central monomer, by extending each arm by 1 ethanol at a time. Extend from
# carb* or branch* Junctions to the matching numbered oxy* Junctions, and then
# finally from a carb* Junction to a "term" Junction.
four_polymer.extend(e1, from_junction_name="branch1", to_junction_name="oxy1")
four_polymer.extend(e2, from_junction_name="branch2", to_junction_name="oxy2")
four_polymer.extend(e3, from_junction_name="branch3", to_junction_name="oxy3")
four_polymer.extend(e4, from_junction_name="branch4", to_junction_name="oxy4")

four_polymer.extend(e1, from_junction_name="carb1", to_junction_name="oxy1")
four_polymer.extend(e2, from_junction_name="carb2", to_junction_name="oxy2")
four_polymer.extend(e3, from_junction_name="carb3", to_junction_name="oxy3")
four_polymer.extend(e4, from_junction_name="carb4", to_junction_name="oxy4")

four_polymer.extend(e1, from_junction_name="carb1", to_junction_name="oxy1")
four_polymer.extend(e2, from_junction_name="carb2", to_junction_name="oxy2")
four_polymer.extend(e3, from_junction_name="carb3", to_junction_name="oxy3")
four_polymer.extend(e4, from_junction_name="carb4", to_junction_name="oxy4")

four_polymer.extend(terminal, from_junction_name="carb1", to_junction_name="term")
four_polymer.extend(terminal, from_junction_name="carb2", to_junction_name="term")
four_polymer.extend(terminal, from_junction_name="carb3", to_junction_name="term")
four_polymer.extend(terminal, from_junction_name="carb4", to_junction_name="term")

# Optionally (but recommended), check the netcharge of the monomers is
# preserved and give it a descriptive name
print(f"netcharge = {four_polymer.topology.netcharge}")
four_polymer.topology.title = "four arm star polymer" # rename your ITP header and image name


# ----- Save the polymer dendrimer to an itp file -----

# Optionally, use Visualize to generate an image of the structure with RDKit
# for an easy visual structure check
four_polymer.save_to_file('data/four_arm_star.json') # optional, text dump in a dictionary format
four_polymer.topology.to_ITP('data/four_arm_star.itp') # write the itp to be used for simulation
Visualize.polymer(four_polymer,infer_bond_order=False).draw2D('data/four_arm_star.png',(400,300)) # optional, visualize the structure in 2D