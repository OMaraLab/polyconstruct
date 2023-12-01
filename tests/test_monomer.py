import json
from pathlib import Path
from polytop.monomer import Monomer
from polytop.topology import Topology
from polytop.junction import Junction, Junctions


def test_monomer(data_dir: Path ):
    arg = Topology.from_ITP(data_dir/"arginine.itp")
    bond_a = arg.get_bond('N3','H20')
    bond_b = arg.get_bond('C11','O1')
    
    monomer = Monomer(arg, [Junction("A",bond_a), Junction("C",bond_b)])
    
    assert len(monomer.topology.atoms) == 26
    assert len(monomer.junctions.get_junctions()) == 2
    assert monomer.junctions.named("A")[0].location.atom_a.atom_name == 'N3'
    assert monomer.junctions.named("A")[0].location.atom_b.atom_name == 'H20'
    assert monomer.junctions.named("C")[0].location.atom_a.atom_name == 'C11'
    assert monomer.junctions.named("C")[0].location.atom_b.atom_name == 'O1'
    
def test_serializable(data_dir: Path, output_dir: Path):
    arg = Topology.from_ITP(data_dir/"arginine.itp")
    bond_a = arg.get_bond('N3','H20')
    bond_b = arg.get_bond('C11','O1')
    junctions = Junctions()
    junctions.add(Junction("A",bond_a))
    junctions.add(Junction("C",bond_b))
 
    
    monomer = Monomer(arg, junctions)
    monomer.save(output_dir/"arg.json")
    
    new_monomer = Monomer.load(output_dir/"arg.json")
    assert len(monomer.topology.atoms) == len(new_monomer.topology.atoms)
    assert monomer.junctions.get_junctions()[0].name == new_monomer.junctions.get_junctions()[0].name
    assert monomer.junctions.get_junctions()[1].name == new_monomer.junctions.get_junctions()[1].name
    assert monomer.junctions.get_junctions()[0].location.atom_a.atom_name == new_monomer.junctions.get_junctions()[0].location.atom_a.atom_name
    assert monomer.junctions.get_junctions()[0].location.atom_b.atom_name == new_monomer.junctions.get_junctions()[0].location.atom_b.atom_name
    assert monomer.junctions.get_junctions()[1].location.atom_a.atom_name == new_monomer.junctions.get_junctions()[1].location.atom_a.atom_name
    assert monomer.junctions.get_junctions()[1].location.atom_b.atom_name == new_monomer.junctions.get_junctions()[1].location.atom_b.atom_name
