from __future__ import annotations

import json
import random
import re
from typing import Dict, List, Optional, Tuple, Union
from polytop.bonds import Bond
from polytop.junction import Junction, Junctions
from polytop.topology import Topology
from polytop.polymer import Polymer
from polytop.atoms import Atom
from polytop.monomer import Monomer
from polytop.visualize import Visualize
import datetime
import copy


class Automatic:
    """
    doctrings go brr
    """
    def __init__(self, orderedList, itp_monomers, numMonomers, monomerList, junctions, dummies) -> None:
        self.directions = orderedList
        self.sizeOfPolymer = numMonomers
        self.monomers = monomerList
        self.junction_atoms = junctions
        self.leaving_atoms = dummies
        self.monomer_paths_dict = itp_monomers

    def build(self, outputName = 'polymer'):
        monomers = []
        names = ["from", "to"]
        print(self.monomer_paths_dict)
        print(self.junction_atoms[0])
        print(self.leaving_atoms[0])
        print(names[0])

        for index, monomer in enumerate(self.directions):
            print(monomer)
            m = Topology.from_ITP(self.monomer_paths_dict[self.directions[index]])
            junctions = []
            if (index == 0): #keeping leaving
                junctions.append(Junction(m.get_atom(self.junction_atoms[0]), m.get_atom(self.leaving_atoms[0]), name = names[0]))
            elif (index == self.sizeOfPolymer-1):
                junctions.append(Junction(m.get_atom(self.junction_atoms[1]), m.get_atom(self.leaving_atoms[1]), name = names[1]))
            else:
                for i, j in enumerate(self.junction_atoms):
                    junctions.append(Junction(m.get_atom(j), m.get_atom(self.leaving_atoms[i]), name = names[i]))
            mono = Monomer(m, junctions)
            monomers.append(mono)
        
        polymer = Polymer(monomers[0])
        for i in range(1, self.sizeOfPolymer):
            polymer.extend(monomers[i], "from", "to")

        polymer.save_to_file(f'{outputName}.json')
        polymer_topology = polymer.topology
        polymer_topology.to_ITP(f'{outputName}.itp')
        Visualize.polymer(polymer,infer_bond_order=False).draw2D(f'{outputName}.png',(400,200))