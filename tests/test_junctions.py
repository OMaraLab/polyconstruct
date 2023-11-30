import json
from polytop.monomer import Monomer
from polytop.topology import Topology
from polytop.junction import Junction, Junctions


def test_junctions():
    arg = Topology.from_ITP("tests/samples/arginine.itp")
    carboxylic_junction = Junction("carboxylic",arg.get_bond('C11','O1'))
    amine1_junction = Junction("amine",arg.get_bond('N3','H20'))
    amine2_junction = Junction("amine",arg.get_bond('N6','H23'))
    junctions = Junctions.from_list([carboxylic_junction, amine1_junction, amine2_junction]) 
    assert len(junctions) == 3
    assert junctions.named("carboxylic") == [carboxylic_junction]
    assert junctions.named("amine") == [amine1_junction, amine2_junction]

def test_junction_serialization():
    arg = Topology.from_ITP("tests/samples/arginine.itp")
    atoms = arg.atoms
    carboxylic_junction = Junction("carboxylic",arg.get_bond('C11','O1'))
    with open("tests/samples/cjunction.json", "w") as f:
        json.dump(carboxylic_junction.to_dict(), f)
    with open("tests/samples/cjunction.json", "r") as f:
        new_junction = Junction.from_dict(json.load(f), atoms)
    assert carboxylic_junction.name == new_junction.name
    assert carboxylic_junction.location.atom_a.atom_name == new_junction.location.atom_a.atom_name
    assert carboxylic_junction.location.atom_b.atom_name == new_junction.location.atom_b.atom_name
    
    
    
def test_junctions_serialization():
    arg = Topology.from_ITP("tests/samples/arginine.itp")
    atoms = arg.atoms
    carboxylic_junction = Junction("carboxylic",arg.get_bond('C11','O1'))
    amine1_junction = Junction("amine",arg.get_bond('N3','H20'))
    amine2_junction = Junction("amine",arg.get_bond('N6','H23'))
    junctions = Junctions()
    junctions.add(carboxylic_junction)
    junctions.add(amine1_junction)
    junctions.add(amine2_junction)
    with open("tests/samples/junctions.json", "w") as f:
        json.dump(junctions.to_dict(), f)
    with open("tests/samples/junctions.json", "r") as f:
        new_junctions = Junctions.from_dict(json.load(f), atoms)
    assert len(new_junctions) == 3
    assert junctions[0].name == new_junctions[0].name
    assert junctions[1].name == new_junctions[1].name
    assert junctions[2].name == new_junctions[2].name
    assert junctions[0].location.atom_a.atom_name == new_junctions[0].location.atom_a.atom_name
    assert junctions[0].location.atom_b.atom_name == new_junctions[0].location.atom_b.atom_name
    assert junctions[1].location.atom_a.atom_name == new_junctions[1].location.atom_a.atom_name
    assert junctions[1].location.atom_b.atom_name == new_junctions[1].location.atom_b.atom_name
    assert junctions[2].location.atom_a.atom_name == new_junctions[2].location.atom_a.atom_name
    assert junctions[2].location.atom_b.atom_name == new_junctions[2].location.atom_b.atom_name
    
