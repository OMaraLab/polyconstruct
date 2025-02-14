# Construction of an ethylamine dendrimer

# import required classes from PolyTop
from polytop.Junction import Junction
from polytop.Monomer import Monomer
from polytop.Visualize import Visualize
from polytop.Polymer import Polymer
from polytop.Topology import Topology

# load in monomer topologies from ITP files
core_mono = Topology.from_ITP("data/dendrimer_core.itp")
bifurcating_mono = Topology.from_ITP("data/dendrimer_bifurcating.itp")
terminal_mono = Topology.from_ITP("data/dendrimer_terminal.itp")

# create junctions for different 'levels' of monomers, depending on what level of branching they exist at in the dendrimer
# junctions are created with the bonding atom and then the leaving atom specified, in that order, and should be given a unique name
central1 = Junction(core_mono.get_atom("N2"), core_mono.get_atom("C9"), name="C1")
central2 = Junction(core_mono.get_atom("N2"), core_mono.get_atom("C7"), name="C2")
central3 = Junction(core_mono.get_atom("N1"), core_mono.get_atom("C2"), name="C3")
central4 = Junction(core_mono.get_atom("N1"), core_mono.get_atom("C3"), name="C4")

b1 = Junction(bifurcating_mono.get_atom("C6"), bifurcating_mono.get_atom("N2"), name = "to")
b2a = Junction(bifurcating_mono.get_atom("N1"), bifurcating_mono.get_atom("C2"), name = "from1")
b2b = Junction(bifurcating_mono.get_atom("N1"), bifurcating_mono.get_atom("C3"), name = "from2")

c1 = Junction(bifurcating_mono.get_atom("C6"), bifurcating_mono.get_atom("N2"), name = "to2")
c2a = Junction(bifurcating_mono.get_atom("N1"), bifurcating_mono.get_atom("C2"), name = "from12")
c2b = Junction(bifurcating_mono.get_atom("N1"), bifurcating_mono.get_atom("C3"), name = "from22")

t = Junction(terminal_mono.get_atom("C1"), terminal_mono.get_atom("N1"), name = "term")

# create monomers from their topologies and any specified junctions
# note that different 'levels' of monomers have different names for their junctions.
# This ensures that layers are added on sequentially. 
# If all monomers have the same junction names the polymerisation will be random and often defaults to a linear shape!
central = Monomer(core_mono, [central1, central2, central3, central4])
bifur1 = Monomer(bifurcating_mono, [b1, b2a, b2b])
bifur2 = Monomer(bifurcating_mono, [c1, c2a, c2b])
cap = Monomer(terminal_mono, [t])

# start the polymer with the central monomer
polymer = Polymer(central)

# extend first layer of dendrimer
polymer.extend(bifur1, from_junction_name="C1", to_junction_name="to")
polymer.extend(bifur1, from_junction_name="C2", to_junction_name="to")
polymer.extend(bifur1, from_junction_name="C3", to_junction_name="to")
polymer.extend(bifur1, from_junction_name="C4", to_junction_name="to")

# extend second layer of dendrimer
# note how the first and second 'layers' have different monomers with different
# Junction names to ensure the correct structure is always formed. Junction and
# Monomer name ambiguity WILL cause formation of random, unreplicatable polymer topologies
polymer.extend(bifur2, from_junction_name="from1", to_junction_name="to2")
polymer.extend(bifur2, from_junction_name="from2", to_junction_name="to2")
polymer.extend(bifur2, from_junction_name="from1", to_junction_name="to2")
polymer.extend(bifur2, from_junction_name="from2", to_junction_name="to2")
polymer.extend(bifur2, from_junction_name="from1", to_junction_name="to2")
polymer.extend(bifur2, from_junction_name="from2", to_junction_name="to2")
polymer.extend(bifur2, from_junction_name="from1", to_junction_name="to2")
polymer.extend(bifur2, from_junction_name="from2", to_junction_name="to2")

# finish polymer by extending on the capping monomer
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

# check netcharge of the monomers is preserved (in this case, close to 0)
print(polymer.topology.netcharge)

# save the dendrimer to a file and visualise the structure with RDKit for an easy visual structure check
polymer.save_to_file('data/ethylamine_dendrimer.json') # text dump
polymer.topology.to_ITP('data/ethylamine_dendrimer.itp')
Visualize.polymer(polymer,infer_bond_order=False).draw2D('data/ethylamine_dendrimer.png',(400,300))