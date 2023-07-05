    
from polytop.monomer import Monomer
from polytop.polymer import Polymer
from polytop.topology import Topology


def test_polymer():    
    arg = Topology.from_ITP("tests/samples/arginine.itp")
    glu = Topology.from_ITP("tests/samples/glutamine.itp")
    bond_a = arg.get_bond('N3','H20')
    bond_b = arg.get_bond('C11','O1')
    arg_monomer = Monomer(arg, bond_a, bond_b)
    
    bond_a = glu.get_bond('N1','H6')
    bond_b = glu.get_bond('C4','O1')
    glu_monomer = Monomer(glu, bond_a, bond_b)
    
    polymer = Polymer([arg_monomer,glu_monomer],['X0','X1'],[20,80], num_monomers= 12, seed= 42, start_monomer= arg_monomer)
    polymer.save_to_file('tests/samples/polymer.json')
    polymer_topology = polymer.get_topology()