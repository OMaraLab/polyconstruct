import json
from polytop.monomer import Monomer
from polytop.topology import Topology
from polytop.polymerization_type import PolymerizationType

def test_PolymerJunction():
    junction_type = PolymerizationType("carboxylic","amino")
    assert junction_type.junction_a == "carboxylic"
    assert junction_type.junction_b == "amino"