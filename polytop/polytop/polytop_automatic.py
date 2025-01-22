from __future__ import annotations

import json
import random
import re
from tqdm import tqdm
from typing import Dict, List, Optional, Tuple, Union
from .Bonds import Bond
from .Junction import Junction, Junctions
from .Topology import Topology
from .Polymer import Polymer
from .Atoms import Atom
from .Monomer import Monomer
from .Visualize import Visualize
import datetime
import copy


class Automatic:
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

        for index, monomer in enumerate(self.directions):
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
        for i in tqdm(range(1, self.sizeOfPolymer), desc='Building polymer topology with PolyTop'):
            polymer.extend(monomers[i], "from", "to")

        polymer.save_to_file(f'{outputName}.json')
        polymer_topology = polymer.topology
        polymer_topology.to_ITP(f'{outputName}.itp')
        try:
            Visualize.polymer(polymer,infer_bond_order=False).draw2D(f'{outputName}.png',(400,200))
        except Exception as e:
            print("Unable to generate RDKit render")
            print(f"Failed with: '{e}'")