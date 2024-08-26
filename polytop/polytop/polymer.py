from __future__ import annotations

import json
import random
import re
from typing import Dict, List, Optional, Tuple, Union
from polytop.bonds import Bond
from polytop.junction import Junction, Junctions
from .topology import Topology
import datetime
import copy

class Polymer:

    def __init__(self, monomer):
        """ 
        A polymer is a topology with a set of junctions that represent the polymerization sites of the monomer.
        Args:
            monomer (Monomer): the monomer to create the polymer from
        """
        new_monomer = copy.deepcopy(monomer)
        self.topology = new_monomer.topology
        self.junctions = new_monomer.junctions

    def has_junction(self, name):
        """Returns True if the polymer has a junction with the given name """
        return any(junction.name == name for junction in self.junctions)
    
    def DFS(self, atom, visited, exclude=None):
        """ 
        Depth first search of the polymer to find all atoms connected to the atom including the atom itself 
        Args:
            atom (Atom): the atom to start the search from
            visited (set): a set of atoms that have already been visited
            exclude (Atom): an atom to exclude from the search (ie: the other side of the junction)
        """
        visited.add(atom)
        for neighbor in atom.bond_neighbours():
            if neighbor is not exclude and neighbor not in visited:
                self.DFS(neighbor, visited, exclude)

    def extend(self, monomer, from_junction_name, to_junction_name, keep_charge = False, bond_length_func = None):
        """
        Extend the polymer by adding a monomer to the polymerization junctions of the polymer
        Args:
            monomer (Monomer): the monomer to add to the polymer
            from_junction_name (str): the name of the junction in the polymer to extend from (a random junction is chosen if multiple with this name are present)
            to_junction_name (str): the name of the junction in the monomer to extend to
            keep_charge (bool): if True, the charge of the polymer is kept the same by forcing the final topology to be 
                the same net charge as the sum charge of the initial topology and the monomer 
            bond_length_func (function): a function that takes 2 bonds (the junction bond on the polymer, and the 
            junction bond on the monomer) and returns a bond length for the junction bond of the extended polymer
        """
        # set up a default bond length function that averages leaving bonds if none is provided
        if bond_length_func is None:
            bond_length_func = lambda bond1, bond2: (bond1.bond_length + bond2.bond_length) / 2


        # if we have keep_charge set then we need to note the existing charges of the polymer and monomer
        monomer_charge = monomer.topology.netcharge
        polymer_charge = self.topology.netcharge
        
        # take a copy of the topology of the monomer 
        new_monomer = monomer.copy()

        # renumber residues in the monomer to avoid conflicts with the polymer
        new_monomer.topology.renumber_residues(self.topology.max_residue_id())
        
        # choose a random polymerization junction of the monomer named to_junction_name to extend this monomer from
        to_junctions = [junction for junction in new_monomer.junctions if junction.name == to_junction_name]
        if not to_junctions:
            raise ValueError(f"No junction named {to_junction_name} found in the monomer")
        to_junction = random.choice(to_junctions)
        if to_junction is None:
            raise ValueError(f"No junction named {to_junction_name} found in the monomer")
        else:
            new_monomer.junctions.remove(to_junction)

        # choose a random polymerization junction of the polymer with from_junction_name to extend this monomer into
        from_junctions = [junction for junction in self.junctions if junction.name == from_junction_name]
        if not from_junctions:
            raise ValueError(f"No junction named {from_junction_name} found in the polymer")
        from_junction = random.choice(from_junctions)
        if from_junction is None:
            raise ValueError(f"No junction named {from_junction_name} found in the polymer")
        else:
            self.junctions.remove(from_junction)

        # calculate the new bond length for the junction bond between the polymer and monomer
        from_junction_bond = self.topology.get_bond(from_junction.monomer_atom, from_junction.residue_atom)
        to_junction_bond = new_monomer.topology.get_bond(to_junction.monomer_atom, to_junction.residue_atom)
        new_bond_length = bond_length_func(from_junction_bond, to_junction_bond)

        # set both junction bonds (polymer and incoming monomer) bond lengths to the new bond length
        from_junction_bond.bond_length = new_bond_length
        to_junction_bond.bond_length = new_bond_length

        # residue atoms other than the first residue atom as identified in the junctions from the polymer and monomer are discarded

        # do a depth first search of the monomer to find all atoms connected to the to_junction 
        discard_from_monomer=set()
        self.DFS(to_junction.residue_atom, discard_from_monomer, exclude=to_junction.monomer_atom)

        # remove the atom in the residue of the junction from the discard_from_monomer set
        discard_from_monomer.remove(to_junction.residue_atom)

        # do a depth first search of the polymer to find all atoms connected to the from_junction including the second atom in the junction
        discard_from_polymer=set()
        self.DFS(from_junction.residue_atom, discard_from_polymer, exclude=from_junction.monomer_atom)

        # remove the atom in the residue of the junction from the discard_from_polymer set
        discard_from_polymer.remove(from_junction.residue_atom)

        for atom in discard_from_polymer:
            self.topology.remove_atom(atom)
        self.topology.reorder_atoms() # renumber the atoms in the polymer so the ids are corrected
        # Note, we still have one atom on both topologies that will be removed to_junction.residue_atom in the monomer and from_junction.residue_atom in the polymer

        # remove the leaving atoms and associated bonds, angles, pairs, and dihedrals from the monomer and polymer
        for atom in discard_from_monomer:
            new_monomer.topology.remove_atom(atom)
        new_monomer.topology.reorder_atoms() # renumber the atoms in the monomer so the ids are corrected
        # renumber all atoms above the max atom id of the polymer setting the formerly attribute of each atom to the old value
        new_monomer.topology.renumber_atoms(max(atom.atom_id for atom in self.topology.atoms))

        # Add the monomer's topology to the polymer
        self.topology.add(new_monomer.topology)

        # Add the monomer's junctions to the polymer fixing up where the id has changed
        for junction in new_monomer.junctions:
            new_monomer_atom = self.topology.get_former_atom(junction.monomer_atom.formerly)
            new_residue_atom = self.topology.get_former_atom(junction.residue_atom.formerly)
            new_junction = Junction(monomer_atom=new_monomer_atom, residue_atom=new_residue_atom, name=junction.name)
            self.junctions.add(new_junction)

        # fix up the to_junction atoms to use the former atom ids
        to_junction.monomer_atom = self.topology.get_former_atom(to_junction.monomer_atom.formerly)
        to_junction.residue_atom = self.topology.get_former_atom(to_junction.residue_atom.formerly)

        # Dihedrals containing a residue atom need to be copied with the remaining atom of the other monomer
        # Angles containing a residue need to be copied with the remaining atom of the other monomer
        # Bonds containing a residue atom need to be copied with the remaining atom of the other monomer
        self.topology.change_atom(from_junction.residue_atom, to_junction.monomer_atom)
        self.topology.change_atom(to_junction.residue_atom, from_junction.monomer_atom)

        # this will double up to_junction.monomer_atom, and from_junction.monomer_atom so we remove the original
        self.topology.atoms.remove(to_junction.monomer_atom) # remove the original to_junction.monomer_atom
        self.topology.atoms.remove(from_junction.monomer_atom) # remove the original from_junction.monomer_atom
        
        # discard both the remaining residue atoms on either side of the junction
        self.topology.remove_atom(from_junction.residue_atom)
        self.topology.remove_atom(to_junction.residue_atom)

        # remove all residue atoms and any bonds angles pair dihedrals containing leaving atoms
        # Any bonds doubling up between 2 atoms need to be averaged and replaced with a single bond

        # deduplicate bonds in the polymer (removing the extra bonds, angles, dihedrals)
        self.topology.deduplicate()

        # reset all the former attributes from the atoms in the polymer
        self.topology.clear_former_ids()

        # if keep_charge is set then we need to adjust the net charge of the polymer to match the sum of the monomer and polymer charges           
        if keep_charge:
            new_charge = monomer_charge + polymer_charge 
            self.topology.netcharge = new_charge

    def save_to_file(self, filename: str) -> None:
        """Save the polymer to a json file"""
        with open(filename, "w") as f:
            json.dump(self.to_dict(), f)

    @classmethod
    def load_from_file(cls, filename: str) -> None:
        """ Load a polymer from a json file"""
        with open(filename, "r") as f:
            data = json.load(f)
        return cls.from_dict(data)

    def to_dict(
        self,
    ) -> Dict[
        str,
        Union[
            List[Dict[str, Union[float, int]]],
            int,
            Optional[Dict[str, Union[float, int]]],
        ],
    ]:
        return {
            "topology": self.topology.to_dict(),
            "junctions": self.junctions.to_dict()
        }

    @classmethod
    def from_dict(
        cls,
        data: Dict[
            str,
            Union[
                List[Dict[str, Union[float, int]]],
                int,
                Optional[Dict[str, Union[float, int]]],
            ],
        ],
    ) -> Polymer:
        new_topology = Topology.from_dict(data["topology"])
        new_junctions = Junctions.from_dict(data["junctions"], new_topology.atoms)
        return cls(
            topology=new_topology,
            junctions=new_junctions,
        )
        
    def __repr__(self) -> str:
        return f"Polymer(({len(self.topology.atoms)} atoms), junctions:{self.junctions.count()})"