import json
from polytop.monomer import Monomer
from polytop.topology import Topology


def test_monomer():
    arg = Topology.from_ITP("tests/samples/arginine.itp")
    monomer = Monomer(arg, arg.get_bond_by_name('N3','H20'), arg.get_bond_by_name('C11','O1'))
    assert len(monomer.LHS.atoms) == 2
    assert len(monomer.link.atoms) == 26-3+2
    assert len(monomer.RHS.atoms) == 3
    
def test_serializable():
    arg = Topology.from_ITP("tests/samples/arginine.itp")
    monomer = Monomer(arg, arg.get_bond(21,20), arg.get_bond(23,25))
    monomer.save("tests/samples/arg.json")
    
    new_monomer = Monomer.load("tests/samples/arg.json")
    assert len(monomer.topology.atoms) == len(new_monomer.topology.atoms)
    assert monomer.bond_a.atom_a.atom_id == new_monomer.bond_a.atom_a.atom_id
    assert monomer.bond_a.atom_b.atom_id == new_monomer.bond_a.atom_b.atom_id
    
    
    
