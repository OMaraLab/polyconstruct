    
from pathlib import Path
import random

import pytest
from polytop.polytop.Junction import Junction
from polytop.polytop.Monomer import Monomer
from polytop.polytop.Visualize import Visualize
from polytop.polytop.Polymer import Polymer
from polytop.polytop.Topology import Topology



def test_simple_polymer(data_dir: Path, output_dir: Path):    
    arg = Topology.from_ITP(data_dir/"arginine.itp")
    arg_N = arg.junction('N3','H20').named("N")
    arg_C = arg.junction('C11','O1').named("C")
    arg_monomer = Monomer(arg, [arg_N, arg_C])
    Visualize.monomer(arg_monomer).draw2D(output_dir/'arginine_monomer.png',(400,200))

    glu = Topology.from_ITP(data_dir/"glutamine.itp")
    glu_N = glu.junction('N1','H6').named("N")
    glu_C = glu.junction('C4','O1').named("C")
    glu_monomer = Monomer(glu, [glu_N, glu_C])
    Visualize.monomer(glu_monomer).draw2D(output_dir/'glutamine_monomer.png',(400,200))

    polymer = Polymer(arg_monomer)

    assert len(arg.atoms) == 26 # arginine has 26 atoms
    assert len(glu.atoms) == 20 # glutamine has 20 atoms
    arg_glu = arg + glu
    assert len(arg_glu.atoms) == 46 # a naive join contains 46 atoms

    polymer.extend(glu_monomer, from_junction_name= 'N', to_junction_name= 'C')
    polymer.topology.title = "arginine glutamine dipeptide"
    polymer.topology.preamble = ["Generated by PolyTop","Author: Richard A. Morris","pytest: tests/test_polymer.py"]    
    
    polymer.save_to_file(output_dir/'simple_polymer.json')
    polymer_topology = polymer.topology
    polymer_topology.to_ITP(output_dir/'simple_polymer.itp')
    Visualize.polymer(polymer).draw2D(output_dir/'simple_polymer.png',(400,200))

    assert polymer.has_junction("N") # we still have an N terminal junction
    assert polymer.has_junction("C") # we still have a C terminal junction
    assert len(polymer.junctions) == 2 # we have 2 junctions left after 2 were joined
    assert len(polymer.topology.atoms) == 43 # we have 43 atoms in the polymer after 3 left
    assert len(polymer.topology.bonds) < len(arg_glu.bonds) # we have fewer bonds than the naive join
    assert len(polymer.topology.angles) < len(arg_glu.angles) # we have fewer angles than the naive join
    
    with pytest.raises(ValueError):
        polymer.topology.get_atom("O2")

    assert polymer.topology.get_atom("O2",1), "we should be able to get the O2 atom from the arg monomer if we specify the residue"
    assert polymer.topology.get_atom("O2",2), "we should be able to get the O2 atom from the glu monomer if we specify the residue"


def test_complex_polymer(data_dir: Path, output_dir: Path):    
    arg = Topology.from_ITP(data_dir/"arginine.itp")
    arg_a1 = arg.junction('N3','H20').named("N")
    arg_a2 = arg.junction('N6','H23').named("N-side-chain")
    arg_c = arg.junction('C11','O1').named("C")

   
    glu = Topology.from_ITP(data_dir/"glutamine.itp")
    glu_N = glu.junction('N1','H6').named("N")
    glu_C = glu.junction('C4','O1').named("C")
    
    arg_monomer = Monomer(arg, [arg_a1, arg_a2, arg_c]) # arginine monomer with 3 junctions
    glu_monomer = Monomer(glu, [glu_N, glu_C]) # glutamine monomer with 2 junctions

    # create a bond length function that returns 10x the average length of the two bonds to the leaving groups
    def long_junction_bonds(bond1, bond2):
        return (bond1.bond_length + bond2.bond_length)*5

    # set random seed for reproducibility
    random.seed(42)
    
    # start with a randomly selected monomer, 20:80::GLU:ARG
    if random.random() < 0.2:
        polymer = Polymer(glu_monomer)
        last_monomer = "glutamine"
    else:
        polymer = Polymer(arg_monomer)
        last_monomer = "arginine"
        
    # add 12 monomers to the polymer backbone alternating         
    for i in range(12):
        if last_monomer == "arginine":
            polymer.extend(glu_monomer, from_junction_name= 'N', to_junction_name= 'C', bond_length_func=long_junction_bonds)
            last_monomer = "glutamine"
        else:
            polymer.extend(arg_monomer, from_junction_name= 'N', to_junction_name= 'C', bond_length_func=long_junction_bonds)
            last_monomer = "arginine"
    
    # extend the backbone in both directions with 2 more arginines
    for i in range(2):
        polymer.extend(arg_monomer, from_junction_name= 'N', to_junction_name= 'C', bond_length_func=long_junction_bonds)
        polymer.extend(arg_monomer, from_junction_name= 'C', to_junction_name= 'N', bond_length_func=long_junction_bonds)
                
    # add another monomer (80:20::GLU:ARG) to the branches of any remaining junction points 3 times
    countdown = 3
    while polymer.has_junction("N-side-chain") and countdown > 0:
        if random.random() < 0.2:
            polymer.extend(arg_monomer, from_junction_name= 'N-side-chain', to_junction_name= 'C', bond_length_func=long_junction_bonds)
        else:
            polymer.extend(glu_monomer, from_junction_name= 'N-side-chain', to_junction_name= 'C', bond_length_func=long_junction_bonds)            
        countdown -= 1
    polymer.topology.title = "complex polymer"
    polymer.topology.preamble = ["Generated by PolyTop","Author: Richard A. Morris","pytest: tests/test_polymer.py"]

    bond_set = set()
    for bond in polymer.topology.bonds:
        bond_set.add(bond)

    assert len(bond_set) == len(polymer.topology.bonds) # no duplicate bonds
        

    polymer.save_to_file(output_dir/'polymer.json')
    polymer_topology = polymer.topology
    polymer_topology.to_ITP(output_dir/'polymer.itp')
    Visualize.polymer(polymer,infer_bond_order=False).draw2D(output_dir/'polymer.png',(400,200))

def test_convert_to_monomer(data_dir:Path):
    glucose_topology = Topology.from_ITP(data_dir/'glucose.itp')

    alpha1 = glucose_topology.junction('C5','O4').named('1')
    alpha6 = glucose_topology.junction('O2','H1').named('6')
    alpha4 = glucose_topology.junction('O3','H2').named('4')
    glucose_topology.title= "Cyclic glucose monomer"

    Glucose_14 = Monomer(glucose_topology, [alpha1, alpha4])

    def create_14_chain(n):
        A_chain = Polymer(Glucose_14)
        for i in range(n-1):
            A_chain.extend(Glucose_14, from_junction_name = "1", to_junction_name = "4")
        return A_chain
    
    def no_duplicate_bonds(t: Topology) -> bool:
        bond_set = set()
        duplicate_bonds = []
        
        for bond in t.bonds:
            bond_tuple = (bond.atom_a.atom_id, bond.atom_b.atom_id)
            if bond_tuple in bond_set:
                duplicate_bonds.append(bond_tuple)
            else:
                bond_set.add(bond_tuple)
        
        if duplicate_bonds:
            print(f"Duplicate bonds found: {duplicate_bonds}")
        
        no_duplicates = len(t.bonds) == len(bond_set)
        return no_duplicates

    chain = create_14_chain(2)
    assert no_duplicate_bonds(chain.topology), "Duplicate bonds in polymer"
    to_monomer = Monomer.from_Polymer(chain)
    assert no_duplicate_bonds(to_monomer.topology), "Duplicate bonds in monomer"


def test_polymer_of_polymers(data_dir: Path, output_dir: Path):
    random.seed(42)
    glucose_topology = Topology.from_ITP(data_dir/'glucose.itp')

    alpha1 = glucose_topology.junction('C5','O4').named('1')
    alpha6 = glucose_topology.junction('O2','H1').named('6')
    alpha4 = glucose_topology.junction('O3','H2').named('4')
    glucose_topology.title= "Cyclic glucose monomer"

    Glucose_14 = Monomer(glucose_topology, [alpha1, alpha4])
    Glucose_146 = Monomer(glucose_topology, [alpha1, alpha4, alpha6])

    def create_14_chain(n):
        A_chain = Polymer(Glucose_14)
        for i in range(n-1):
            A_chain.extend(Glucose_14, from_junction_name = "1", to_junction_name = "4")
        return A_chain

    def create_146_chain(n):
        B_chain = Polymer(Glucose_146)
        for i in range(n-1):
            B_chain.extend(Glucose_146, from_junction_name = "1", to_junction_name = "4")
        return B_chain



    # Construct 2 B chains with 4 A chains branched from each B chain
    B_chains=[]
    for i in range(2):
        chain = create_146_chain(9)
        for j in range(4):
            A_chain = create_14_chain(9)
            A_chain_as_monomer = Monomer.from_Polymer(A_chain)
            chain.extend(A_chain_as_monomer, from_junction_name = "6", to_junction_name = "1")
        chain.topology.title = f"Beta chain {i}"
        chain_as_monomer = Monomer.from_Polymer(chain)
        B_chains.append(chain_as_monomer)

    C_chain = create_146_chain(8)
    C_chain.topology.title = "Glycogen"
    Glycogen = C_chain
    for monomer in B_chains:
        Glycogen.extend(monomer, from_junction_name = "6", to_junction_name = "1")

    Visualize.polymer(Glycogen).draw2D(output_dir/'glycogen.png',(400,200))
    
def test_multijoined_polymer(data_dir: Path, output_dir: Path):   
    glu_ltop = Topology.from_ITP(data_dir/"glucose.itp")
    glu_lj1 = Junction(glu_ltop.get_atom("C6"), glu_ltop.get_atom("O2"), name = "O2")
    glu_lj2 = Junction(glu_ltop.get_atom("O3"), glu_ltop.get_atom("H2"), name = "O3")
    glu_left = Monomer(glu_ltop, [glu_lj1, glu_lj2])

    glu_rtop = Topology.from_ITP(data_dir/"glucose.itp")
    glu_rj1 = Junction(glu_rtop.get_atom("O5"), glu_rtop.get_atom("H4"), name = "O5")
    glu_rj2 = Junction(glu_rtop.get_atom("C4"), glu_rtop.get_atom("O6"), name = "O6")
    glu_right = Monomer(glu_rtop, [glu_rj1, glu_rj2])
    
    polymer = Polymer(glu_left)
    polymer.extend(glu_right, from_junction_name="O2", to_junction_name="O5")
    polymer.extra_bond(from_junction_name="O3", to_junction_name="O6")
                
    polymer.topology.title = "multijoined glucose"
    polymer.topology.preamble = ["Generated by PolyTop","Author: Richard A. Morris","pytest: tests/test_polymer.py"]

    bond_set = set()
    for bond in polymer.topology.bonds:
        bond_set.add(bond)

    assert len(bond_set) == len(polymer.topology.bonds) # no duplicate bonds

    assert len(polymer.topology.bonds) == 44
        
    polymer.save_to_file(output_dir/'double_sugar_polymer.json')
    polymer_topology = polymer.topology
    polymer_topology.to_ITP(output_dir/'double_sugar_polymer.itp')
    Visualize.polymer(polymer,infer_bond_order=False).draw2D(output_dir/'double_sugar_polymer.png',(400,200))

def test_OPLS_polymer(data_dir: Path, output_dir: Path):
    """
    This test is to make sure a Polymer constructed from OPLS ff monomers
    extends without error
    """
    top = Topology.from_ITP(data_dir/"OPLS_UNK_460A12.itp", format="opls")

    # Create a Junction to join 'to' and another to join 'from'.
    # Provide the bonding atom and the leaving atom, in that order, for the
    # Junction - they must have a bond between them.
    to_j = Junction(top.get_atom("C02"), top.get_atom("C01"), name = "to")
    from_j = Junction(top.get_atom("C03"), top.get_atom("C04"), name = "from")

    # Create a Monomer from the Topology and a list of the Junctions
    monomer = Monomer(top, [to_j, from_j])

    # Start the Polymer with one Monomer
    polymer = Polymer(monomer)

    # Extend the Polymer to the desired length (in this case 20)
    for i in range(29):
        polymer.extend(monomer, from_junction_name="from", to_junction_name="to")

    # Save the polymer to a file and visualise the structure with RDKit for an easy visual structure check
    polymer.topology.to_ITP(output_dir/'OPLS_30mer.itp')
    Visualize.polymer(polymer,infer_bond_order=False).draw2D(output_dir/'OPLS_30mer.png',(400,300))

def test_AMBER_polymer(data_dir: Path, output_dir: Path):
    """
    This test is to make sure a Polymer constructed from AMBER ff monomers
    extends without error
    """
    top = Topology.from_ITP(data_dir/"AMBER_PNIPAM_extend.itp", format="amber")

    # Create a Junction to join 'to' and another to join 'from'.
    # Provide the bonding atom and the leaving atom, in that order, for the
    # Junction - they must have a bond between them.
    to_j = Junction(top.get_atom("C02"), top.get_atom("C01"), name = "to")
    from_j = Junction(top.get_atom("C03"), top.get_atom("C04"), name = "from")

    # Create a Monomer from the Topology and a list of the Junctions
    monomer = Monomer(top, [to_j, from_j])

    # Start the Polymer with one Monomer
    polymer = Polymer(monomer)

    # Extend the Polymer to the desired length (in this case 20)
    for i in range(29):
        polymer.extend(monomer, from_junction_name="from", to_junction_name="to")

    # Save the polymer to a file and visualise the structure with RDKit for an easy visual structure check
    polymer.topology.to_ITP(output_dir/'AMBER_30mer.itp')
    Visualize.polymer(polymer,infer_bond_order=False).draw2D(output_dir/'AMBER_30mer.png',(400,300))

def test_CHARMM_polymer(data_dir: Path, output_dir: Path):
    """
    This test is to make sure a Polymer constructed from CHARMM ff monomers
    extends without error
    """
    top = Topology.from_ITP(data_dir/"CHARMM_PNIPAM_extend.itp", format="charmm")

    # Create a Junction to join 'to' and another to join 'from'.
    # Provide the bonding atom and the leaving atom, in that order, for the
    # Junction - they must have a bond between them.
    to_j = Junction(top.get_atom("C2"), top.get_atom("C1"), name = "to")
    from_j = Junction(top.get_atom("C3"), top.get_atom("C4"), name = "from")

    # Create a Monomer from the Topology and a list of the Junctions
    monomer = Monomer(top, [to_j, from_j])

    # Start the Polymer with one Monomer
    polymer = Polymer(monomer)

    # Extend the Polymer to the desired length (in this case 20)
    for i in range(29):
        polymer.extend(monomer, from_junction_name="from", to_junction_name="to")

    # Save the polymer to a file and visualise the structure with RDKit for an easy visual structure check
    polymer.topology.to_ITP(output_dir/'CHARMM_30mer.itp')
    Visualize.polymer(polymer,infer_bond_order=False).draw2D(output_dir/'CHARMM_30mer.png',(400,300))