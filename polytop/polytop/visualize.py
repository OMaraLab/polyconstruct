import copy
import logging
import re
from typing import Tuple
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw
from polytop.junction import Junctions
from polytop.topology import Topology
from functools import singledispatchmethod
from polytop.polymer import Polymer
from polytop.monomer import Monomer

class Visualize:
    def __init__(
        self,
        topology: Topology,
        junctions: Junctions = None,
        infer_bond_order = True,
    ):
        self.topology = copy.deepcopy(topology)
        self.junctions = Junctions() if junctions is None else junctions
        self.atom_mapping = {}

        if infer_bond_order:
            self.infer_bond_order()
        
    def infer_bond_order(self):
        # add double bonds to the topology using the known valencies of the atoms and number of declared bonds
        self.valencies = {
            "H": 1,
            "F": 1,
            "Cl": 1,  # this can also be 3,5,7 depending on the oxidation state of elements it's bonded to
            "Br": 1,  # this can also be 3,5,7 depending on the oxidation state of elements it's bonded to
            "I": 1,  # this can also be 3,5,7 depending on the oxidation state of elements it's bonded to
            "O": 2,
            "S": 2,
            "N": 3,  # this can also be 2,4,5 depending on the oxidation state of elements it's bonded to
            "C": 4,
            "X": 1,
        }

        # Iterate over the atoms and count the number of bonds for each atom
        bond_count = {atom: 0 for atom in self.topology.atoms}
        for bond in self.topology.bonds:
            bond_count[bond.atom_a] += 1
            bond_count[bond.atom_b] += 1

        # Compare the bond count with the known valency and infer bond orders
        inferred_bond_orders = []
        for bond in self.topology.bonds:
            atom_a = bond.atom_a
            atom_b = bond.atom_b
            element_a = atom_a.element
            element_b = atom_b.element

            if element_a in self.valencies and element_b in self.valencies:
                valency_a = self.valencies[element_a]
                valency_b = self.valencies[element_b]

                bond_count_a = bond_count[atom_a]
                bond_count_b = bond_count[atom_b]

                if bond_count_a < valency_a and bond_count_b < valency_b:
                    inferred_bond_order = 2
                else:
                    inferred_bond_order = 1

                inferred_bond_orders.append((atom_a, atom_b, inferred_bond_order))
                   
        # Update the bond orders in the Topology object
        for atom_a, atom_b, bond_order in inferred_bond_orders:
            self.topology.get_bond(atom_a.atom_id, atom_b.atom_id).order = bond_order
                    
            
    @classmethod
    def polymer(cls, polymer: Polymer, infer_bond_order = True):
        return cls(topology=polymer.topology, junctions=polymer.junctions, infer_bond_order = infer_bond_order)

    @classmethod
    def monomer(cls, monomer: Monomer, infer_bond_order = True):
        return cls(topology=monomer.topology, junctions=monomer.junctions, infer_bond_order = infer_bond_order)

    @classmethod
    def topology(cls, topology: Topology, infer_bond_order = True):
        return cls(topology=topology, infer_bond_order = infer_bond_order)

    def to_rdKit_Chem_mol(self,options: dict[str,any]):

        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)
        logging.basicConfig(level=logging.ERROR)
        mol = Chem.RWMol()

        for atom in self.topology.atoms:
            atom_name = atom.atom_name
            # if the atom is a virtual atom (X), then set the element to H for rdKit 
            element = atom.element
            new_atom = Chem.Atom(element)
            self.atom_mapping[atom.atom_id] = mol.AddAtom(new_atom)

        for bond in self.topology.bonds:
            try:
                if bond.order == 2:
                    bond_type = Chem.BondType.DOUBLE
                elif bond.order == 3:
                    bond_type = Chem.BondType.TRIPLE
                elif bond.order == 1:
                    bond_type = Chem.BondType.SINGLE
                else:
                    bond_type = Chem.BondType.UNSPECIFIED

                mol.AddBond(
                    self.atom_mapping[bond.atom_a.atom_id],
                    self.atom_mapping[bond.atom_b.atom_id],
                    bond_type,
                )
            except KeyError as e:
                print(f"KeyError for bond: {bond}, KeyError: {e}")

        # Sanitize the molecule without virtual atoms
        Chem.SanitizeMol(mol)

        # Label atoms
        for atom in self.topology.atoms:
            mol_atom = mol.GetAtomWithIdx(self.atom_mapping[atom.atom_id])
            atom_name = atom.atom_name
            element = atom.element
            index = atom.index
            atom_label=''
            if options["show_atom_ID"]:
                atom_label = f"{element}<sub>{index}</sub>"
            else:
                if not element == "C": # don't label carbons
                    atom_label = element
                    if options["remove_explicit_H"]:
                        count_h = sum([1 for neighbour in atom.bond_neighbours() if neighbour.element == "H"])
                        if count_h == 1:
                            atom_label += f"H"
                        elif count_h > 1:
                            atom_label += f"H<sub>{count_h}</sub>"
            mol_atom.SetProp("atomLabel", atom_label)
        
        if options.get("highlight_junctions",False):
            # Draw bonds that represent junction locations in a different colour
            for junction in self.junctions:
                atom_a = junction.monomer_atom
                atom_b = junction.residue_atom
                bond = mol.GetBondBetweenAtoms(self.atom_mapping[atom_a.atom_id], self.atom_mapping[atom_b.atom_id])
                index = bond.GetIdx()
                bond.SetProp("Junction", junction.name)
                bond.SetProp("bondNote", '"'+junction.name+'"')
                atom = mol.GetAtomWithIdx(self.atom_mapping[atom_b.atom_id])
                atom.SetProp("Junction", junction.name)

        return mol

    def draw3D(self, view=None, highlight_junctions: bool = False,show_atom_ID:bool= False, remove_explicit_H:bool = True):
        options = {
            "highlight_junctions":highlight_junctions, 
            "show_atom_ID":show_atom_ID,
            "remove_explicit_H": remove_explicit_H
            }
        mol = self.to_rdKit_Chem_mol(options=options)  # Removed the MolToMolBlock conversion
        Chem.SanitizeMol(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        view.addModel(Chem.MolToMolBlock(mol), "mol")  # Convert to MolBlock for py3Dmol
        view.setStyle({"stick": {}})
        view.zoomTo()
        
    def remove_non_junction_hydrogens(self, mol):
        # Get the indices of all hydrogen atoms
        hydrogen_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1]

        # Get the indices of the hydrogens connected to a bond with the "Junction" property
        junction_hydrogen_indices = [bond.GetBeginAtomIdx() for bond in mol.GetBonds() if bond.HasProp("Junction") and bond.GetBeginAtom().GetAtomicNum() == 1]
        junction_hydrogen_indices.extend(bond.GetEndAtomIdx() for bond in mol.GetBonds() if bond.HasProp("Junction") and bond.GetEndAtom().GetAtomicNum() == 1)

        # Remove the indices of the junction hydrogens from the list of hydrogen indices
        hydrogen_indices = [idx for idx in hydrogen_indices if idx not in junction_hydrogen_indices]

        # Sort the indices in descending order
        hydrogen_indices.sort(reverse=True)

        # Remove the hydrogens
        mol = Chem.RWMol(mol)
        for idx in hydrogen_indices:
            mol.RemoveAtom(idx)

        return mol.GetMol()

    def draw2D(
        self,
        filename: str,
        size: Tuple[int,int] = (600,300),
        highlight_junctions = False,
        show_atom_ID = False,
        show_legend = True,
        remove_explicit_H = True,
    ):
        """
        Draw the molecule using RDKit's DrawMolecule function
        Args:
            filename (str): the filename to save the image to
            size (tuple[int,int]): the size of the image to draw
            highlight_junctions (bool): if True, highlight the junctions in the molecule
            show_atom_ID (bool): if True, show the atom ID in the molecule
            show_legend (bool): if True, show the legend in the molecule
            remove_explicit_H (bool): if True, remove the explicit hydrogens from the molecule        
        """
        options = {
            "size":size, 
            "highlight_junctions":highlight_junctions, 
            "show_atom_ID":show_atom_ID,
            "show_legend":show_legend,
            "remove_explicit_H":remove_explicit_H,
            }
        if options["show_atom_ID"] and options["remove_explicit_H"]:
            raise ValueError("show_atom_ID and remove_explicit_H cannot both be True")  
        
        mol = self.to_rdKit_Chem_mol(options)
        
        # Initialize args and kwargs for DrawMolecule
        draw_kwargs = {}
        
        if options.get("remove_explicit_H",True):
            mol = self.remove_non_junction_hydrogens(mol)
        
        size = options.get("size",(600,300))
        d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1]) 
        opt = d.drawOptions()

        if options.get("highlight_junctions",False):
            junction_color = (0,0,1,.2)
            text_color = (0,0,0,1)
            opt.minFontSize = 16
            opt.highlightBondWidthMultiplier = 16
            opt.annotationFontScale = 1.2
            opt.setAnnotationColour(text_color)
            opt.setHighlightColour(junction_color)
            
            mol_bonds = mol.GetBonds()
            junction_bonds = [bond.GetIdx() for bond in mol_bonds if bond.HasProp("Junction")]
            junction_bonds_colors = {bond: junction_color for bond in junction_bonds}

            junction_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.HasProp("Junction")]
            junction_atoms = []
            draw_kwargs["highlightAtoms"] = junction_atoms
            draw_kwargs["highlightBonds"] = junction_bonds
            draw_kwargs["highlightBondColors"] = junction_bonds_colors
            draw_kwargs["highlightAtomColors"] = None

        d.SetFontSize(0.75)
        # Draw the molecule
        if options.get("show_legend",True):
            rdMolDraw2D.PrepareAndDrawMolecule(d, mol, legend=self.topology.title, **draw_kwargs)
        else:
            d.DrawMolecule(mol, **draw_kwargs)
    
        # Finish drawing
        d.FinishDrawing()

        # Save the image
        with open(filename, 'wb') as f:
            f.write(d.GetDrawingText())


    