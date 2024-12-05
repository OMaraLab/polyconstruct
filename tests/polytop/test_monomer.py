import json
from pathlib import Path
from polytop.monomer import Monomer
from polytop.topology import Topology
from polytop.junction import Junction, Junctions


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

def test_renumber_monomer(data_dir: Path):
    glu = Topology.from_ITP(data_dir/"glucose.itp")
    glu_1 = glu.junction("C1","O3").named("1")
    glu_4 = glu.junction("C4","HC3").named("4")
    monomer = Monomer(glu, [glu_1, glu_4])

    monomer.renumber_atoms(100)
    assert monomer.topology.atoms[0].atom_id == 101
    assert monomer.topology.atoms[1].atom_id == 102
    assert len(monomer.junctions) == 2
    assert monomer.junctions.get_junctions()[0].monomer_atom.atom_name == 'C1'
    assert monomer.junctions.get_junctions()[0].residue_atom.atom_name == 'O3'
    assert monomer.junctions.get_junctions()[1].monomer_atom.atom_name == 'C4'
    assert monomer.junctions.get_junctions()[1].residue_atom.atom_name == 'HC3'

