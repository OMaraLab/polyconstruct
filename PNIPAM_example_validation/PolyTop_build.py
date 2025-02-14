from polytop.Junction import Junction
from polytop.Monomer import Monomer
from polytop.Polymer import Polymer
from polytop.Topology import Topology

top = Topology.from_ITP("monomers/UNK_460A12_PolyTop.itp", format="opls")

# Create a Junction to join 'to' and another to join 'from'.
# Provide the bonding atom and the leaving atom, in that order, for the
# Junction - they must have a bond between them.
from_j = Junction(top.get_atom("C02"), top.get_atom("C01"), name = "from")
to_j = Junction(top.get_atom("C03"), top.get_atom("C04"), name = "to")

# Create a Monomer from the Topology and a list of the Junctions
monomer = Monomer(top, [to_j, from_j])

# Start the Polymer with one Monomer
polymer = Polymer(monomer)

# Extend the Polymer to the desired length (in this case 30)
for i in range(29):
    polymer.extend(monomer, from_junction_name="to", to_junction_name="from")

# Save the polymer to a file and visualise the structure with RDKit for an easy visual structure check
polymer.topology.title = "UNK_460A12" # renames the ITP header and image
polymer.topology.to_ITP('polymers/30mer.itp')