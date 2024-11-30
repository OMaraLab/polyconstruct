import json
from pathlib import Path
from polytop.monomer import Monomer
from polytop.topology import Topology
from polytop.junction import Junction, Junctions


def test_junctions(data_dir: Path):
    arg = Topology.from_ITP(data_dir / "arginine.itp")
    # 3 different ways of creating a junction
    amine1_junction = Junction(arg.get_atom('N6'), arg.get_atom('H24'),"amine")
    amine2_junction = arg.junction('N6','H23').named("amine")
    amine3_junction = Junction(arg.get_atom('N5'), arg.get_atom('H26'))
    carboxylic_junction = Junction.from_topology(arg,'C11','O1').named("carboxylic")
    junctions = Junctions([carboxylic_junction, amine1_junction, amine2_junction, amine3_junction])
    assert len(junctions) == 4
    assert junctions.named("carboxylic") == [carboxylic_junction]
    assert junctions.named("amine") == [amine1_junction, amine2_junction]
    assert junctions.named("N5-H26") == [amine3_junction]
    assert junctions[0].monomer_atom.atom_name == 'C11'
    assert junctions[0].residue_atom.atom_name == 'O1'
    assert junctions[1].monomer_atom.atom_name == 'N6'
    assert junctions[1].residue_atom.atom_name == 'H24'
    assert junctions[2].monomer_atom.atom_name == 'N6'
    assert junctions[2].residue_atom.atom_name == 'H23'
    assert junctions[3].monomer_atom.atom_name == 'N5'
    assert junctions[3].residue_atom.atom_name == 'H26'

def test_junction_serialization(data_dir: Path, output_dir: Path):
    arg = Topology.from_ITP(data_dir / "arginine.itp")
    atoms = arg.atoms
    carboxylic_junction = arg.junction('C11','O1').named("carboxylic")
    
    # Write to the file
    (output_dir / "junction.json").write_text(json.dumps(carboxylic_junction.to_dict()))

    # Read from the file
    new_junction = Junction.from_dict(json.loads((output_dir / "junction.json").read_text()),atoms)
    
    assert carboxylic_junction.name == new_junction.name
    assert carboxylic_junction.monomer_atom.atom_name == new_junction.monomer_atom.atom_name
    assert carboxylic_junction.residue_atom.atom_name == new_junction.residue_atom.atom_name
    
    
    
def test_junctions_serialization(data_dir: Path, output_dir: Path):
    arg = Topology.from_ITP(data_dir/"arginine.itp")
    atoms = arg.atoms
    carboxylic_junction = arg.junction('C11','O1').named("carboxylic")
    amine1_junction = arg.junction("N6","H23").named("amine")
    amine2_junction = arg.junction("N6","H24").named("amine")
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
    assert junctions[0].monomer_atom.atom_name == new_junctions[0].monomer_atom.atom_name
    assert junctions[1].monomer_atom.atom_name == new_junctions[1].monomer_atom.atom_name
    assert junctions[2].monomer_atom.atom_name == new_junctions[2].monomer_atom.atom_name
    assert junctions[0].residue_atom.atom_name == new_junctions[0].residue_atom.atom_name
    assert junctions[1].residue_atom.atom_name == new_junctions[1].residue_atom.atom_name
    
    
    
