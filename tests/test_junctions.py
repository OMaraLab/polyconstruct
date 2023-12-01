import json
from pathlib import Path
from polytop.monomer import Monomer
from polytop.topology import Topology
from polytop.junction import Junction, Junctions


def test_junctions(data_dir: Path):
    arg = Topology.from_ITP(data_dir / "arginine.itp")
    carboxylic_junction = Junction("carboxylic",arg.get_bond('C11','O1'))
    amine1_junction = Junction("amine",arg.get_bond('N3','H20'))
    amine2_junction = Junction("amine",arg.get_bond('N6','H23'))
    junctions = Junctions.from_list([carboxylic_junction, amine1_junction, amine2_junction]) 
    assert len(junctions) == 3
    assert junctions.named("carboxylic") == [carboxylic_junction]
    assert junctions.named("amine") == [amine1_junction, amine2_junction]

def test_junction_serialization(data_dir: Path, output_dir: Path):
    arg = Topology.from_ITP(data_dir / "arginine.itp")
    atoms = arg.atoms
    carboxylic_junction = Junction("carboxylic",arg.get_bond('C11','O1'))
    
    # Write to the file
    (output_dir / "junction.json").write_text(json.dumps(carboxylic_junction.to_dict()))

    # Read from the file
    new_junction = Junction.from_dict(json.loads((output_dir / "junction.json").read_text()),atoms)
    
    assert carboxylic_junction.name == new_junction.name
    assert carboxylic_junction.location.atom_a.atom_name == new_junction.location.atom_a.atom_name
    assert carboxylic_junction.location.atom_b.atom_name == new_junction.location.atom_b.atom_name
    
    
    
def test_junctions_serialization(data_dir: Path, output_dir: Path):
    arg = Topology.from_ITP(data_dir/"arginine.itp")
    atoms = arg.atoms
    carboxylic_junction = Junction("carboxylic",arg.get_bond('C11','O1'))
    amine1_junction = Junction("amine",arg.get_bond('N3','H20'))
    amine2_junction = Junction("amine",arg.get_bond('N6','H23'))
    junctions = Junctions()
    junctions.add(carboxylic_junction)
    junctions.add(amine1_junction)
    junctions.add(amine2_junction)

    # Write to the file
    (output_dir / "junctions.json").write_text(json.dumps(junctions.to_dict()))

    # Read from the file
    new_junctions = Junctions.from_dict(json.loads((output_dir / "junctions.json").read_text()),atoms)

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
    
