from __future__ import annotations

import json
import random
from typing import Dict, List, Optional, Tuple, Union, Callable
from .pairs import Pair
from .bonds import Bond
from .atoms import Atom
from .junction import Junction, Junctions
from .dihedrals import Dihedral
from .angles import Angle
from .topology import Topology
import copy

# degenerate functions that follow these prototypes are passed to Polymer.extend so the caller can override the default behaviors.
# these functions are passed the topology, the from and to polymerization junctions and a list of the degenerate bonds
# angles, or dihedrals that are about to be removed
# The function should add back any bonds, angles, or dihedrals to the topology that are required to maintain the structure of the polymer

DegenerateBondFunc = Callable[[Bond, Bond, Topology, Atom, Atom ], None] 
DegenerateAngleFunc = Callable[[List[Angle], Topology,Junction, Junction], None]
DegenerateDihedralFunc = Callable[[List[Dihedral], Topology,Junction, Junction], None] 

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
        self.joined_monomer = None

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

    def remove_all_but_2_Leaving_atoms(self, from_junction, to_junction, adding_monomer = True):
        # leaving atoms other than the first and second leaving atom as identified in the junctions from the polymer and monomer are discarded

        # do a depth first search of the monomer to find all atoms connected to the to_junction 
        discard_from_monomer=set()
        self.DFS(to_junction.leaving_atom, discard_from_monomer, exclude=to_junction.remaining_atom)

        # remove the leaving atom of the junction from the discard_from_monomer set
        discard_from_monomer.remove(to_junction.leaving_atom)
        if to_junction.second_leaving_atom is not None:
            discard_from_monomer.remove(to_junction.second_leaving_atom)

        # do a depth first search of the polymer to find all atoms connected to the from_junction including the second atom in the junction
        discard_from_polymer=set()
        self.DFS(from_junction.leaving_atom, discard_from_polymer, exclude=from_junction.remaining_atom)

        # remove the atom in the residue of the junction from the discard_from_polymer set
        discard_from_polymer.remove(from_junction.leaving_atom)
        if from_junction.second_leaving_atom is not None:
            discard_from_polymer.remove(from_junction.second_leaving_atom)

        for atom in discard_from_polymer:
            self.topology.remove_atom(atom)
        self.topology.reorder_atoms() # renumber the atoms in the polymer so the ids are corrected
        # Note, we still have two atoms on both topologies that will be removed to_junction.leaving_atom and to_junction.second_leaving_atom in the monomer 
        # and from_junction.residue_atom and from_junction.second_residue_atom in the polymer

        # remove the leaving atoms and associated bonds, angles, pairs, and dihedrals from the monomer and polymer
        if adding_monomer:
            for atom in discard_from_monomer:
                self.joined_monomer.topology.remove_atom(atom)
            self.joined_monomer.topology.reorder_atoms() # renumber the atoms in the monomer so the ids are corrected
            # renumber all atoms above the max atom id of the polymer setting the formerly attribute of each atom to the old value
            self.joined_monomer.topology.renumber_atoms(max(atom.atom_id for atom in self.topology.atoms))
        else:
            for atom in discard_from_monomer:
                self.topology.remove_atom(atom)
            self.topology.reorder_atoms() # renumber the atoms in the monomer so the ids are corrected
    
    
    def extend(self, monomer, from_junction_name, to_junction_name, keep_charge = True, 
               bond_func:DegenerateBondFunc = None, angle_func:DegenerateAngleFunc = None, dihedral_func:DegenerateDihedralFunc = None):
        """
        Extend the polymer by adding a monomer to the polymerization junctions of the polymer
        Args:
            monomer (Monomer): the monomer to add to the polymer
            from_junction_name (str): the name of the junction in the polymer to extend from (a random junction is chosen if multiple with this name are present)
            to_junction_name (str): the name of the junction in the monomer to extend to
            keep_charge (bool): default True - the charge of the polymer is kept the same by forcing the final topology to be 
                the same net charge as the sum charge of the initial topology and the monomer 
            bond_func (function): default None - a function to manually add the new bond after degenerate bonds are removed
            angle_func (function): default None - a function to manually add new angles after degenerate angles are removed
            dihedral_func (function): default None - a function to manually add new dihedrals after degenerate dihedrals are removed
        """
        def default_bond_func(polymer_bond: Bond, monomer_bond: Bond, topology: Topology, polymer_atom : "Atom", monomer_atom : "Atom"):
            ''' A default bond function that averages the bond length and force constant of the last 2 leaving bonds '''
            new_bond = Bond(polymer_atom, monomer_atom)
            new_bond.bond_type = polymer_bond.bond_type
            new_bond.force_constant = (polymer_bond.force_constant + monomer_bond.force_constant) / 2
            new_bond.bond_length = (polymer_bond.bond_length + monomer_bond.bond_length) / 2

        if bond_func is None:
            bond_func = default_bond_func

        def default_angle_func(polymer_angles: List[Angle], topology: Topology, from_junction: Junction, to_junction: Junction):
            ''' A default angle function that changes any angle that contains the leaving atoms to the new atoms '''
            for angle in polymer_angles:
                angle.replace_atom(from_junction.leaving_atom, to_junction.remaining_atom)
                angle.replace_atom(to_junction.leaving_atom, from_junction.remaining_atom)

        if angle_func is None:
            angle_func = default_angle_func

        def default_dihedral_func(polymer_dihedrals: List[Dihedral], topology: Topology, from_junction: Junction, to_junction: Junction):
            for dihedral in polymer_dihedrals:
                dihedral.replace_atom(from_junction.leaving_atom, to_junction.remaining_atom)
                dihedral.replace_atom(to_junction.leaving_atom, from_junction.remaining_atom)
                dihedral.replace_atom(from_junction.second_leaving_atom, to_junction.second_remaining_atom)
                dihedral.replace_atom(to_junction.second_leaving_atom, from_junction.second_remaining_atom)

        if dihedral_func is None:
            dihedral_func = default_dihedral_func

        # if we have keep_charge set then we need to note the existing charges of the polymer and monomer
        monomer_charge = monomer.topology.netcharge
        polymer_charge = self.topology.netcharge
        
        # take a copy of the topology of the monomer 
        self.joined_monomer = monomer.copy()

        # renumber residues in the monomer to avoid conflicts with the polymer
        self.joined_monomer.topology.renumber_residues(self.topology.max_residue_id())
        
        # choose a random polymerization junction of the monomer named to_junction_name to extend this monomer from
        to_junctions = [junction for junction in self.joined_monomer.junctions if junction.name == to_junction_name]
        if not to_junctions:
            raise ValueError(f"No junction named {to_junction_name} found in the monomer")
        to_junction = random.choice(to_junctions)
        if to_junction is None:
            raise ValueError(f"No junction named {to_junction_name} found in the monomer")
        else:
            self.joined_monomer.junctions.remove(to_junction)

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
        from_junction_bond = self.topology.get_bond(from_junction.remaining_atom, from_junction.leaving_atom)
        to_junction_bond = self.joined_monomer.topology.get_bond(to_junction.remaining_atom, to_junction.leaving_atom)
        bond_func(from_junction_bond, to_junction_bond, self.topology, from_junction.remaining_atom, to_junction.remaining_atom)

        # Discard leaving atoms and reorder atoms
        self.remove_all_but_2_Leaving_atoms(from_junction, to_junction)


        # Add the monomer's topology to the polymer
        self.topology.add(self.joined_monomer.topology)

        # calculate the new angles for the junction angles between the polymer and monomer
        angles = []
        for angle in self.topology.angles:
            if angle.references_atom(from_junction.leaving_atom):
                angles.append(angle) 
            elif angle.references_atom(to_junction.leaving_atom):
                angles.append(angle)
        angle_func(angles, self.topology, from_junction, to_junction)

        dihedrals = []
        for dihedral in self.topology.dihedrals:
            if dihedral.references_atom(from_junction.leaving_atom):
                dihedrals.append(dihedral) 
            elif dihedral.references_atom(to_junction.leaving_atom):
                dihedrals.append(dihedral)
            elif dihedral.references_atom(from_junction.second_leaving_atom):
                dihedrals.append(dihedral)
            elif dihedral.references_atom(to_junction.second_leaving_atom):
                dihedrals.append(dihedral)

        dihedral_func(dihedrals, self.topology, from_junction, to_junction)

        # Add the monomer's junctions to the polymer fixing up where the id has changed
        for junction in self.joined_monomer.junctions:
            new_monomer_atom = self.topology.get_former_atom(junction.monomer_atom.formerly)
            new_residue_atom = self.topology.get_former_atom(junction.residue_atom.formerly)
            new_junction = Junction(remaining_atom=new_monomer_atom, leaving_atom=new_residue_atom, name=junction.name)
            self.junctions.add(new_junction)

        # fix up the to_junction atoms to use the former atom ids
        to_junction.monomer_atom = self.topology.get_former_atom(to_junction.monomer_atom.formerly)
        to_junction.residue_atom = self.topology.get_former_atom(to_junction.residue_atom.formerly)


        # discard both the remaining residue atoms on either side of the junction
        if from_junction.second_remaining_atom is not None:
            self.topology.remove_atom(from_junction.second_remaining_atom)
        if to_junction.second_residue_atom is not None:
            self.topology.remove_atom(to_junction.second_residue_atom)
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
        return f"Polymer(({len(self.topology.atoms)} atoms), junctions:{len(self.junctions)})"