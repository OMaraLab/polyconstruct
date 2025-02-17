#!/usr/bin/env python 
import pytest
import random
from polyconf.polyconf.Polymer import Polymer
from polyconf.polyconf.Monomer import Monomer
import MDAnalysis as mda
from .test_data import test_pdbs, PEI_pdbs, PEI_start, PEI_monomer, PEI_end

@pytest.mark.parametrize("pdb_file", test_pdbs)
def test_polymer_initialization(data_dir, pdb_file):
    '''
    Test that we can initialize a Polymer object from a pdb file
    And that it defines a bunch of attributes and functions that we expect
    '''
    
    first_monomer = Monomer(f"{data_dir}/{pdb_file}")
    polymer = Polymer(first_monomer)
    assert polymer.first == first_monomer
    assert polymer.first.residues.resids[0] == 1

    # let's make sure some attributes we expect are set
    assert hasattr(polymer, 'first')
    assert hasattr(polymer, 'polymer')
    assert hasattr(polymer, 'atoms')

    # let's make sure some functions we expect are defined
    methods = ['select_atoms', 'copy', 'renamer', 'extend', 'dist', 'gencomp',
               'shuffle', 'shuffler', 'dihedral_solver', 'gen_pairlist',
               'newresid', 'maxresid']
    for method in methods:
        assert hasattr(Polymer, method), f"Polymer class should have a method named {method}"

def test_polymer_initialization_arguments():
    '''
    Test that the Polymer class raises an error when initialized with a non-Monomer object
    '''
    with pytest.raises(ValueError):
        Polymer(None)
    with pytest.raises(TypeError):
        Polymer("not a Monomer")

@pytest.mark.parametrize("pdb_file", test_pdbs)
def test_polymer_attributes(data_dir, pdb_file):
    first_monomer = Monomer(f"{data_dir}/{pdb_file}")
    polymer = Polymer(first_monomer)
    assert polymer.first.residues.resids[0] == 1

@pytest.mark.parametrize("pdb_file", test_pdbs)
def test_select_atoms(data_dir, pdb_file: str):
    """
    Test that we can find an atom by name in the polymer using the select_atoms method
    """
    polymer = Polymer(Monomer(data_dir / pdb_file))
    random_atom = random.choice(polymer.atoms)
    selected_atoms = polymer.select_atoms(f"name {random_atom.name}")
    assert len(selected_atoms) == 1, f"There should only be one {random_atom.name} atom in {pdb_file}"
    assert selected_atoms[0].name == random_atom.name, f"The selected atom should be named {random_atom.name}"

@pytest.mark.parametrize("pdb_file", test_pdbs)
def test_polymer_copy(data_dir, pdb_file: str):
    """
    Test that the Polymer copy() function creates a new, identical Polymer.
    """
    polymer = Polymer(Monomer(data_dir / pdb_file))
    new_polymer = Polymer.copy(polymer)
    assert len(polymer.select_atoms("all")) == len(new_polymer.select_atoms("all"))
    for atom in polymer.atoms:
        assert len(new_polymer.select_atoms(f"name {atom.name}")) == 1
    assert polymer.polymer.dimensions == new_polymer.polymer.dimensions
    random_atom = random.choice(polymer.atoms)
    selected_atom = polymer.select_atoms(f"name {random_atom.name}")[0]
    selected_atom.name = "different"
    assert len(polymer.select_atoms(f"name different")) == 1
    assert len(new_polymer.select_atoms(f"name different")) == 0 # changing atom in original should not change atom in the copy

@pytest.mark.parametrize("pdb_file", test_pdbs)
def test_polymer_renamer(data_dir, pdb_file: str):
    """
    Test that the Polymer renamer() function renames the correct atoms to the
    provided new name.
    """
    polymer = Polymer(Monomer(data_dir / pdb_file))
    name = random.choice(polymer.atoms).name
    polymer.renamer(1, name, nameout='X')
    assert len(polymer.select_atoms(f"name X*")) == 1
    name2 = random.choice(polymer.atoms).name
    polymer.renamer(1, name2, nameout='DUMMY')
    assert len(polymer.select_atoms(f"name DUMMY1")) == 1

def test_polymer_extend(data_dir):
    """
    Test that the Polymer extend works as expected.
    """
    polymer = Polymer(Monomer(data_dir / 'PEI_start.pdb'))

    adds=128
    for i in range (0,adds):
        polymer.extend( # extend with one monomer, aligned along this step's linearization vector
                Monomer(data_dir / 'PEI_monomer.pdb'), # extend with this monomer
                n=polymer.maxresid(), # we will allways add onto the existing monomer with the highest resid
                nn=polymer.newresid(), # the incoming monomer needs a new resid
                names=dict(Q='CX',P='C1',S='N1',R='NX'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
                joins=[('N1','C1')],# new connection between N1_i and C1_i+1 
                ) 

    polymer.extend(
                Monomer(data_dir / 'PEI_end.pdb'),
                n=polymer.maxresid(), 
                nn=polymer.newresid(),
                names=dict(Q='CX',P='C1',S='N1',R='NX'),
                joins=[('N1','C1')],
                )
    
    assert polymer.maxresid() == 130

def test_polymer_maxresid(data_dir):
    """
    Test that the Polymer maxresid() function returns the correct number.
    """
    polymer = Polymer(Monomer(data_dir / 'PEI_start.pdb'))
    assert polymer.maxresid() == 1

    adds=8
    for i in range (0,adds):
        polymer.extend( # extend with one monomer, aligned along this step's linearization vector
                Monomer(data_dir / 'PEI_monomer.pdb'), # extend with this monomer
                n=polymer.maxresid(), # we will allways add onto the existing monomer with the highest resid
                nn=polymer.newresid(), # the incoming monomer needs a new resid
                names=dict(Q='CX',P='C1',S='N1',R='NX'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
                joins=[('N1','C1')],# new connection between N1_i and C1_i+1 
                ) 

    polymer.extend(
                Monomer(data_dir / 'PEI_end.pdb'),
                n=polymer.maxresid(), 
                nn=polymer.newresid(),
                names=dict(Q='CX',P='C1',S='N1',R='NX'),
                joins=[('N1','C1')],
                )
    
    assert polymer.maxresid() == 10

    new_polymer = Polymer(Monomer.monomer_from_u(polymer.polymer), keep_resids=False)
    assert new_polymer.maxresid() == 1 # if keep_resids==False, the new polymer will only have resid 1

def test_polymer_newresid(data_dir):
    """
    Test that the Polymer newresid() function returns the correct number.
    """
    polymer = Polymer(Monomer(data_dir / 'PEI_start.pdb'))
    assert polymer.newresid() == 2

    adds=8
    for i in range (0,adds):
        polymer.extend( # extend with one monomer, aligned along this step's linearization vector
                Monomer(data_dir / 'PEI_monomer.pdb'), # extend with this monomer
                n=polymer.maxresid(), # we will allways add onto the existing monomer with the highest resid
                nn=polymer.newresid(), # the incoming monomer needs a new resid
                names=dict(Q='CX',P='C1',S='N1',R='NX'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
                joins=[('N1','C1')],# new connection between N1_i and C1_i+1 
                )
        assert polymer.newresid() == i+3

    polymer.extend(
                Monomer(data_dir / 'PEI_end.pdb'),
                n=polymer.maxresid(), 
                nn=polymer.newresid(),
                names=dict(Q='CX',P='C1',S='N1',R='NX'),
                joins=[('N1','C1')],
                )
    
    assert polymer.newresid() == 11

def test_polymer_dihedralsolver_and_genpairlist(data_dir):
    """
    Test that the Polymer dihedral_solver() and gen_pairlist() work as expected
    and are able to resolve the polymer dihedrals.
    """
    polymer = Polymer(Monomer(data_dir / 'PEI_start.pdb'))
    adds=128
    for i in range (0,adds):
        polymer.extend( # extend with one monomer, aligned along this step's linearization vector
                Monomer(data_dir / 'PEI_monomer.pdb'), # extend with this monomer
                n=polymer.maxresid(), # we will allways add onto the existing monomer with the highest resid
                nn=polymer.newresid(), # the incoming monomer needs a new resid
                names=dict(Q='CX',P='C1',S='N1',R='NX'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
                joins=[('N1','C1')],# new connection between N1_i and C1_i+1 
                ) 

    polymer.extend(
                Monomer(data_dir / 'PEI_end.pdb'),
                n=polymer.maxresid(), 
                nn=polymer.newresid(),
                names=dict(Q='CX',P='C1',S='N1',R='NX'),
                joins=[('N1','C1')],
                )

    CN_dihedrals=polymer.gen_pairlist(J='C2',K='N1',first_resid=1,last_resid=130,mult=3)

    cutoff = 0.7

    polymer.dihedral_solver(CN_dihedrals,dummies='CX,NX',cutoff=cutoff)

    for dh in CN_dihedrals:
        # check that the dihedral_solver has removed all clashes
        assert polymer.dist('C2',dh['J_resid'],'N1',dh['K_resid'],'CX,NX',backwards_only=False) >= cutoff

def test_polymer_shuffler(data_dir):
    """
    Test that the Polymer shuffler() works as expected.
    """
    polymer = Polymer(Monomer(data_dir / 'PEI_start.pdb'))
    adds=128
    for i in range (0,adds):
        polymer.extend( # extend with one monomer, aligned along this step's linearization vector
                Monomer(data_dir / 'PEI_monomer.pdb'), # extend with this monomer
                n=polymer.maxresid(), # we will allways add onto the existing monomer with the highest resid
                nn=polymer.newresid(), # the incoming monomer needs a new resid
                names=dict(Q='CX',P='C1',S='N1',R='NX'), # C1_i+1 fit to CX_i, then rotate so NX_i+1 fit to N1_i 
                joins=[('N1','C1')],# new connection between N1_i and C1_i+1 
                ) 

    polymer.extend(
                Monomer(data_dir / 'PEI_end.pdb'),
                n=polymer.maxresid(), 
                nn=polymer.newresid(),
                names=dict(Q='CX',P='C1',S='N1',R='NX'),
                joins=[('N1','C1')],
                )

    CN_dihedrals=polymer.gen_pairlist(J='C2',K='N1',first_resid=1,last_resid=130,mult=3)
    
    distances_before = []

    for dh in CN_dihedrals:
        distances_before.append(polymer.dist('C2',dh['J_resid'],'N1',dh['K_resid'],'CX,NX',backwards_only=False))
    
    polymer.shuffler(CN_dihedrals)

    for i, dh in enumerate(CN_dihedrals):
        # check that the shuffler has chaged all the dihedrals
        assert polymer.dist('C2',dh['J_resid'],'N1',dh['K_resid'],'CX,NX',backwards_only=False) != distances_before[i]

# TODO: Add more tests for other methods and functionalities of the Polymer class as needed