import json
from pathlib import Path
from polytop.monomer import Monomer
from polytop.topology import Topology
from polytop.junction import Junction, Junctions
from polytop.visualize import Visualize


def test_monomer_ARG(data_dir: Path ):
    ARG = Topology.from_ITP(data_dir/"arginine.itp")
    ARG_N = ARG.junction('N3','H20').named("N")
    ARG_C = ARG.junction('C11','O1').named("C")
    ARG_monomer = Monomer(ARG, [ARG_N, ARG_C])
    
    assert len(ARG_monomer.topology.atoms) == 26
    assert len(ARG_monomer.junctions.get_junctions()) == 2
    assert ARG_monomer.junctions.named("N")[0].monomer_atom.atom_name == 'N3'
    assert ARG_monomer.junctions.named("N")[0].residue_atom.atom_name == 'H20'
    assert ARG_monomer.junctions.named("C")[0].monomer_atom.atom_name == 'C11'
    assert ARG_monomer.junctions.named("C")[0].residue_atom.atom_name == 'O1'

def test_monomer_GLN(data_dir: Path ):
    GLN = Topology.from_ITP(data_dir/"glutamine.itp")
    GLN_N = GLN.junction('N1','H6').named("N")
    GLN_C = GLN.junction('C4','O1').named("C")
    
    GLN_monomer = Monomer(GLN, [GLN_N, GLN_C])
    
    assert len(GLN_monomer.topology.atoms) == 20
    assert len(GLN_monomer.junctions.get_junctions()) == 2
    assert GLN_monomer.junctions.named("N")[0].monomer_atom.atom_name == 'N1'
    assert GLN_monomer.junctions.named("N")[0].residue_atom.atom_name == 'H6'
    assert GLN_monomer.junctions.named("C")[0].monomer_atom.atom_name == 'C4'
    assert GLN_monomer.junctions.named("C")[0].residue_atom.atom_name == 'O1'
    
def test_monomer_GLU(data_dir: Path ):
    glu = Topology.from_ITP(data_dir/"glucose.itp")
    
    glu_1 = glu.junction("C1","O3").named("1")
    glu_4 = glu.junction("C4","HC3").named("4")
    glu_6 = glu.junction("C6","HC11").named("6")
    monomer = Monomer(glu, [glu_1, glu_4, glu_6])
    
    assert len(monomer.topology.atoms) == 24
    assert len(monomer.junctions.get_junctions()) == 3
    assert monomer.junctions.named("1")[0].monomer_atom.atom_name == 'C1'
    assert monomer.junctions.named("1")[0].residue_atom.atom_name == 'O3'
    assert monomer.junctions.named("4")[0].monomer_atom.atom_name == 'C4'
    assert monomer.junctions.named("4")[0].residue_atom.atom_name == 'HC3'
    assert monomer.junctions.named("6")[0].monomer_atom.atom_name == 'C6'
    assert monomer.junctions.named("6")[0].residue_atom.atom_name == 'HC11'
    
def test_serializable(data_dir: Path, output_dir: Path):
    arg = Topology.from_ITP(data_dir/"arginine.itp")
    junctions = Junctions()
    junctions.add(arg.junction('N3','H20').named("A"))
    junctions.add(arg.junction('C11','O1').named("C"))
    
    monomer = Monomer(arg, junctions)
    monomer.save(output_dir/"arg.json")
    
    new_monomer = Monomer.load(output_dir/"arg.json")
    assert len(monomer.topology.atoms) == len(new_monomer.topology.atoms)
    assert monomer.junctions.get_junctions()[0].name == new_monomer.junctions.get_junctions()[0].name
    assert monomer.junctions.get_junctions()[1].name == new_monomer.junctions.get_junctions()[1].name
    assert monomer.junctions.get_junctions()[0].monomer_atom.atom_name == new_monomer.junctions.get_junctions()[0].monomer_atom.atom_name
    assert monomer.junctions.get_junctions()[0].residue_atom.atom_name == new_monomer.junctions.get_junctions()[0].residue_atom.atom_name
    assert monomer.junctions.get_junctions()[1].monomer_atom.atom_name == new_monomer.junctions.get_junctions()[1].monomer_atom.atom_name
    assert monomer.junctions.get_junctions()[1].residue_atom.atom_name == new_monomer.junctions.get_junctions()[1].residue_atom.atom_name

def test_renumber_monomer_and_change_residue_id(data_dir: Path, output_dir: Path):
    glu = Topology.from_ITP(data_dir/"glutamine.itp")
    glu_N = glu.junction('N1','H6').named("N")
    glu_C = glu.junction('C4','O1').named("C")
    glu_monomer = Monomer(glu, [glu_N, glu_C])
    Visualize.monomer(glu_monomer).draw2D(output_dir/'glutamine_before_renumber.png',(400,200))
    glu_copy = glu_monomer.copy()
    glu_copy.topology.renumber_atoms(100)
    glu_copy.set_residue_id(100)
    glu_copy.topology.to_ITP(output_dir/'glutamine_renumbered.itp')
    Visualize.monomer(glu_copy).draw2D(output_dir/'glutamine_after_renumber.png',(400,200))    
    for i,atom in enumerate(glu_copy.topology.atoms):
        assert atom.atom_id == glu_monomer.topology.atoms[i].atom_id + 100
        assert atom.residue_id == 100
        assert atom.residue_name == glu_monomer.topology.atoms[i].residue_name
        assert atom.atom_name == glu_monomer.topology.atoms[i].atom_name
        assert atom.atom_type == glu_monomer.topology.atoms[i].atom_type
        assert atom.charge_group_num == glu_monomer.topology.atoms[i].charge_group_num
        assert atom.partial_charge == glu_monomer.topology.atoms[i].partial_charge
        assert atom.mass == glu_monomer.topology.atoms[i].mass
        for bond in atom.bonds:
            assert glu_monomer.topology.get_bond(bond.atom_a.atom_id-100, bond.atom_b.atom_id-100)
        
