# Construction of a 4-arm PEG star polymer from single monomeric units

# import required classes from PolyTop
from polytop.Junction import Junction
from polytop.Monomer import Monomer
from polytop.Visualize import Visualize
from polytop.Polymer import Polymer
from polytop.Topology import Topology

# load in monomer topologies from ITP files
ethanol = Topology.from_ITP("data/extended_ethanol.itp") # main arm monomer
methane = Topology.from_ITP("data/extended_methane.itp") # terminal monomer
neopentane = Topology.from_ITP("data/extended_neopentane.itp") # central monomer

# create junctions for each monomer with the bonding atom and then the leaving atom specified, in that order, with a unique name
oxy_j1 = Junction(ethanol.get_atom("O1"), ethanol.get_atom("C1"), name = "oxy1")
carb_j1 = Junction(ethanol.get_atom("C3"), ethanol.get_atom("O2"), name = "carb1")
oxy_j2 = Junction(ethanol.get_atom("O1"), ethanol.get_atom("C1"), name = "oxy2")
carb_j2 = Junction(ethanol.get_atom("C3"), ethanol.get_atom("O2"), name = "carb2")
oxy_j3 = Junction(ethanol.get_atom("O1"), ethanol.get_atom("C1"), name = "oxy3")
carb_j3 = Junction(ethanol.get_atom("C3"), ethanol.get_atom("O2"), name = "carb3")
oxy_j4 = Junction(ethanol.get_atom("O1"), ethanol.get_atom("C1"), name = "oxy4")
carb_j4 = Junction(ethanol.get_atom("C3"), ethanol.get_atom("O2"), name = "carb4")

j1 = Junction(neopentane.get_atom("C1"), neopentane.get_atom("O1"), name = "branch1")
j2 = Junction(neopentane.get_atom("C3"), neopentane.get_atom("O2"), name = "branch2")
j3 = Junction(neopentane.get_atom("C4"), neopentane.get_atom("O3"), name = "branch3")
j4 = Junction(neopentane.get_atom("C5"), neopentane.get_atom("O4"), name = "branch4")

term_j = Junction(methane.get_atom("C1"), methane.get_atom("O1"), name = "term")

# create monomers from their topologies and any specified junctions
e1 = Monomer(ethanol, [oxy_j1, carb_j1])
e2 = Monomer(ethanol, [oxy_j2, carb_j2])
e3 = Monomer(ethanol, [oxy_j3, carb_j3])
e4 = Monomer(ethanol, [oxy_j4, carb_j4])

central = Monomer(neopentane, [j1, j2, j3, j4])

terminal = Monomer(methane, [term_j]) # only needs one junction to join to the ends of each arm

# start the polymer with the central monomer
four_polymer = Polymer(central)

# attach three ethanols to each of the four junctions (j1-j4) of the central monomer
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

# check polymer charge and give it a descriptive name
print(f"netcharge = {four_polymer.topology.netcharge}")
four_polymer.topology.title = "four arm star polymer - overlapped monomers" # rename your ITP header and image name

# save the polymer to a file and visualise the structure with RDKit for an easy visual structure check
four_polymer.save_to_file('data/four_arm_star_overlapped_monomers.json') # text dump
four_polymer.topology.to_ITP('data/four_arm_star_overlapped_monomers.itp')
Visualize.polymer(four_polymer,infer_bond_order=False).draw2D('data/four_arm_star_overlapped_monomers.png',(400,300))