# Construction of a simple linear homopolymer of PEI

# Import required classes from PolyTop
from polytop.Junction import Junction
from polytop.Monomer import Monomer
from polytop.Visualize import Visualize
from polytop.Polymer import Polymer
from polytop.Topology import Topology

# Load in monomer Topology from ITP file
top = Topology.from_ITP("data/pei.itp")

# Create a Junction to join 'to' and another to join 'from'.
# Provide the bonding atom and the leaving atom, in that order, for the
# Junction - they must have a bond between them.
to_j = Junction(top.get_atom("C51"), top.get_atom("C62"), name = "to")
from_j = Junction(top.get_atom("N7"), top.get_atom("C6"), name = "from")

# Create a Monomer from the Topology and a list of the Junctions
monomer = Monomer(top, [to_j, from_j])

# Start the Polymer with one Monomer
polymer = Polymer(monomer)

# Extend the Polymer to the desired length (in this case 20)
for i in range(19):
    polymer.extend(monomer, from_junction_name="from", to_junction_name="to")

# Save the polymer to a file and visualise the structure with RDKit for an easy visual structure check
polymer.topology.title = "pei polymer" # renames the ITP header and image
polymer.save_to_file('data/pei_linear_polymer.json') # text dump
polymer.topology.to_ITP('data/pei_linear_polymer.itp')
Visualize.polymer(polymer,infer_bond_order=False).draw2D('data/pei_linear_polymer.png',(400,300))