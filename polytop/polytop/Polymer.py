from __future__ import annotations

import json
import random
import re
from typing import Dict, List, Optional, Tuple, Union
from .Bonds import Bond
from .Junction import Junction, Junctions
from .Topology import Topology
import datetime
import copy

class Polymer:
    """
    Represents a polymer, which is a Topology with a set of Junctions, in a
    molecular system.

    :param monomer: the Monomer to create the Polymer from.
    :type monomer: Monomer
    """
    def __init__(self, monomer):
        """
        Represents a polymer, which is a Topology with a set of Junctions, in a
        molecular system.

        :param monomer: the Monomer to create the Polymer from.
        :type monomer: Monomer
        """
        new_monomer = copy.deepcopy(monomer)
        self.topology = new_monomer.topology
        self.junctions = new_monomer.junctions
        self.joined_monomer = None

    def has_junction(self, name: str) -> bool:
        """
        Check if the Polymer contains a Junction with the given name. 

        :param name: the unique name of a Junction.
        :type name: str
        :return: True if the Polymer has a junction with the given name, or
                False if not.
        :rtype: bool
        """
        return any(junction.name == name for junction in self.junctions)
    
    def DFS(self, atom, visited: set, exclude=None):
        """
        Depth first search of the Polymer to find all Atoms connected to the
        Atom provided in 'atom', including the Atom itself. 

        :param atom: the Atom to start the search from.
        :type atom: Atom
        :param visited: a set of atoms that have already been visited, and are
                connected to the Atom 'atom'.
        :type visited: set
        :param exclude: an Atom to exclude from the search (i.e. the other side
                of the junction), defaults to None.
        :type exclude: Atom, optional
        """
        visited.add(atom)
        for neighbor in atom.bond_neighbours():
            if neighbor is not exclude and neighbor not in visited:
                self.DFS(neighbor, visited, exclude)

    def _join_here(self, from_junction, to_junction, adding_monomer = True):
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
    
    def extra_bond(self, from_junction_name: str, to_junction_name: str, bond_length_func = None):
        """
        Form an extra bond within the existing Polymer WITHOUT adding a new
        Monomer. Useful for creating double bonds or complex Polymer
        architectures.

        :param from_junction_name: unique name of a Junction in the Polymer to
                create a new Bond from.
        :type from_junction_name: str
        :param to_junction_name: unique name of a Junction in the Polymer to
                create a new Bond to.
        :type to_junction_name: str
        :param bond_length_func: a lambda function to calculate the length of
                the new Bond from the Bonds in the two provided Junctions,
                defaults to None (such that the new Bond length is the average
                of the two existing Junction's Bonds)
        :type bond_length_func: a lambda function, optional
        :raises ValueError: if either from_junction_name or to_junction_name do
                not correspond to a Junction present in the Polymer with the
                same name.
        """
        # set up a default bond length function that averages leaving bonds if none is provided
        if bond_length_func is None:
            bond_length_func = lambda bond1, bond2: (bond1.bond_length + bond2.bond_length) / 2
        
        # choose a random polymerization junction of the polymer named to_junction_name to create new bond to
        to_junctions = [junction for junction in self.junctions if junction.name == to_junction_name]
        if not to_junctions:
            raise ValueError(f"No junction named {to_junction_name} found in the polymer")
        to_junction = random.choice(to_junctions)
        if to_junction is None:
            raise ValueError(f"No junction named {to_junction_name} found in the polymer")
        else:
            self.junctions.remove(to_junction)

        # choose a random polymerization junction of the polymer with from_junction_name to creat new bond from
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
        to_junction_bond = self.topology.get_bond(to_junction.monomer_atom, to_junction.residue_atom)
        new_bond_length = bond_length_func(from_junction_bond, to_junction_bond)

        # set both junction bonds (polymer and incoming monomer) bond lengths to the new bond length
        from_junction_bond.bond_length = new_bond_length
        to_junction_bond.bond_length = new_bond_length

        # Discard leaving atoms and reorder atoms
        self._join_here(from_junction, to_junction, adding_monomer=False)

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


    def extend(self, monomer, from_junction_name: str, to_junction_name: str, keep_charge = True, bond_length_func = None):
        """
        Extend the Polymer by joining a new Monomer to a Junction of the
        Polymer, via formation of a new Bond between the incoming Monomer's
        'to_junction' and the Polymer's 'from_junction' monomer atoms.

        When using this function, ensure that the Monomer and Polymer Junctions
        being connected have a matching two-atom depth overlap (i.e. use
        extended topologies). This will ensure that the new Bond, Angles and
        Dihedrals are correct.
        
        For best results:
        * Provide topologies where all atoms that will be lost have a partial
        charge of 0. This will ensure the charge of each atom in the resulting
        Polymer is identical to it's charge in the provided Monomer it came from.
        * Ensure that all Junctions within the Polymer, and by extension the
        incoming Monomer, have unique names (unless they are redundant). If
        more than one Junction is present with a given name, one of these
        Junctions will be selected at random for the Polymer extension. This
        will result in unexpected Polymer topologies and prevent reproducibility.

        :param monomer: the Monomer to add to the Polymer
        :type monomer: Monomer
        :param from_junction_name: unique name of a Junction in the Polymer to
                create a new Bond from (i.e. extend from).
        :type from_junction_name: str
        :param to_junction_name: unique name of a Junction in the Monomer to
                create a new Bond to (i.e. extend to).
        :type to_junction_name: str
        :param keep_charge: if True, maintains the net charge of the Polymer by
                redistributing the charges of any lost atoms among the Atoms in
                the two most recently added residues of the Polymer (i.e.
                typically the incoming Monomer and the monomeric residue it is
                joining to), defaults to True
        :type keep_charge: bool, optional
        :param bond_length_func: a lambda function that takes 2 Bonds (the
                Polymer's Junction Bond and the Monomer's Junction Bond) and
                returns a Bond length for the new Bond, defaults to None (such
                that the new Bond length is the average of the two existing
                Junction's Bonds)
        :type bond_length_func: a lambda function, optional
        :raises ValueError: if from_junction_name does not correspond to a
                Junction present in the Polymer with the same name.
        :raises ValueError: if to_junction_name does not correspond to a
                Junction present in the Monomer with the same name.
        """
        # set up a default bond length function that averages leaving bonds if none is provided
        if bond_length_func is None:
            bond_length_func = lambda bond1, bond2: (bond1.bond_length + bond2.bond_length) / 2


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
        from_junction_bond = self.topology.get_bond(from_junction.monomer_atom, from_junction.residue_atom)
        to_junction_bond = self.joined_monomer.topology.get_bond(to_junction.monomer_atom, to_junction.residue_atom)
        new_bond_length = bond_length_func(from_junction_bond, to_junction_bond)

        # set both junction bonds (polymer and incoming monomer) bond lengths to the new bond length
        from_junction_bond.bond_length = new_bond_length
        to_junction_bond.bond_length = new_bond_length

        # Discard leaving atoms and reorder atoms
        self._join_here(from_junction, to_junction)

        # Add the monomer's topology to the polymer
        self.topology.add(self.joined_monomer.topology)

        # Add the monomer's junctions to the polymer fixing up where the id has changed
        for junction in self.joined_monomer.junctions:
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
        """
        Save and export the Polymer to a JSON text dump.

        :param file_path: path to and desired name of output file.
        :type file_path: str
        """
        with open(filename, "w") as f:
            json.dump(self.to_dict(), f)

    @classmethod
    def load_from_file(cls, filename: str) -> Polymer:
        """
        Load a JSON text dump of a Polymer, generated with
        Polymer.save_to_file(), to a new Polymer.

        :param filename: path to the Polymer JSON dictionary text dump file.
        :type filename: str
        :return: a new Polymer
        :rtype: Polymer
        """
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
        """
        Convert this Polymer to a dictionary representation.

        The structure of the dictionary is as below:
        {"topology": self.topology.to_dict(),
        "junctions": self.junctions.to_dict()}

        :return: a dictionary containing references to the dictionary
                representations of this Polymer's Topology and Junctions
                attributes.
        :rtype: Dict[ str, Union[ List[Dict[str, Union[float, int]]], int, Optional[Dict[str, Union[float, int]]], ], ]
        """
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
        """
        Create a new Polymer from a dictionary, such as that created with
        Polymer.to_dict().

        The structure of the dictionary is as below:
        {"topology": self.topology.to_dict(),
        "junctions": self.junctions.to_dict()}

        :param data: dictionary containing data to make a Polymer, generate
                with 'to_dict()'.
        :type data: Dict[ str, Union[ List[Dict[str, Union[float, int]]], int, Optional[Dict[str, Union[float, int]]], ], ]
        :return: a new Polymer
        :rtype: Polymer
        """
        new_topology = Topology.from_dict(data["topology"])
        new_junctions = Junctions.from_dict(data["junctions"], new_topology.atoms)
        return cls(
            topology=new_topology,
            junctions=new_junctions,
        )
        
    def __repr__(self) -> str:
        return f"Polymer(({len(self.topology.atoms)} atoms), junctions:{len(self.junctions)})"