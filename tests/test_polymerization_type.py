import json
from polytop.Monomer import Monomer
from polytop.Topology import Topology
from polytop.Polymerization_type import PolymerizationType

def test_PolymerJunction():
    junction_type = PolymerizationType("carboxylic","amino")
    assert junction_type.junction_a == "carboxylic"
    assert junction_type.junction_b == "amino"