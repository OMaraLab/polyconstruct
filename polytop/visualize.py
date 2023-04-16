import logging
import re
from typing import Tuple

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw

from polytop.topology import Topology


class Visualize:
    def __init__(
        self,
        topology: Topology,
    ):
        self.topology = topology
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
        }

    def infer_bond_orders(self) -> "Visualize":
        # Step 1: Get a dictionary of known valencies for common atom types

        # Step 2: Iterate over the atoms and count the number of bonds for each atom
        bond_count = {atom: 0 for atom in self.topology.atoms}
        for bond in self.topology.bonds:
            bond_count[bond.atom_a] += 1
            bond_count[bond.atom_b] += 1

        # Step 3: Compare the bond count with the known valency and infer bond orders
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

        # Step 4: Update the bond orders in the Topology object
        for atom_a, atom_b, bond_order in inferred_bond_orders:
            self.topology.get_bond(atom_a.atom_id, atom_b.atom_id).order = bond_order
        return self

    def to_rdKit_Chem_mol(self):

        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)
        logging.basicConfig(level=logging.ERROR)
        mol = Chem.RWMol()
        atom_mapping = {}


        for atom in self.topology.atoms:
            atom_name = atom.atom_name
            # Extract the element symbol from the atom type
            element = re.sub("[^a-zA-Z]", "", atom_name)
            new_atom = Chem.Atom(element)
            atom_mapping[atom.atom_id] = mol.AddAtom(new_atom)

        for bond in self.topology.bonds:
            try:
                if bond.order == 2:
                    bond_type = Chem.rdchem.BondType.DOUBLE
                elif bond.order == 3:
                    bond_type = Chem.rdchem.BondType.TRIPLE
                else:
                    bond_type = Chem.rdchem.BondType.SINGLE
                mol.AddBond(
                    atom_mapping[bond.atom_a.atom_id],
                    atom_mapping[bond.atom_b.atom_id],
                    bond_type,
                )
            except KeyError as e:
                print(f"KeyError for bond: {bond}, KeyError: {e}")
                
        for atom in self.topology.pseudoatoms:
            pseudoatom = mol.GetAtomWithIdx(atom_mapping[atom.atom_id])
            pseudoatom.SetAtomicNum(0)
            pseudoatom.SetProp("RDKIT_ATOM_SYMBOL", "X")
            pseudoatom.SetFormalCharge(0)

        Chem.SanitizeMol(mol)

            

        return mol

    def create_py3Dmol_view(self, view=None, show_hydrogens=False):
        # render display presentation of glutamine
        mol = self.to_rdKit_Chem_mol()  # Removed the MolToMolBlock conversion
        Chem.SanitizeMol(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        view.addModel(Chem.MolToMolBlock(mol), "mol")  # Convert to MolBlock for py3Dmol
        view.setStyle({"stick": {}})
        view.zoomTo()

    def create_2D_image(
        self,
        filename: str,
        size: Tuple[int, int] = (600, 300),
        remove_explicit_Hs: bool = True,
    ):

        mol = self.to_rdKit_Chem_mol()
        if remove_explicit_Hs:
            mol = Chem.RemoveHs(mol)
        Chem.SanitizeMol(mol)

        for atom in mol.GetAtoms():
            # Check the number of heavy atoms bonded to the current atom
            heavy_count = 0
            neighbor_symbols = []
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() != "H":
                    heavy_count += 1
                    neighbor_symbols.append(neighbor.GetSymbol())

            # Modify the label of the current atom based on the number of heavy neighbors
            if atom.GetSymbol() == "C":
                atom.SetProp("atomLabel", "")
            elif heavy_count > 1:
                atom.SetProp("atomLabel", atom.GetSymbol())
            else:
                num_h = atom.GetTotalNumHs()
                if num_h == 0:
                    atom.SetProp("atomLabel", atom.GetSymbol())
                elif num_h == 1:
                    atom.SetProp("atomLabel", f"{atom.GetSymbol()}H")
                else:
                    atom.SetProp(
                        "atomLabel",
                        f"{atom.GetSymbol()}H<sub>{atom.GetTotalNumHs()}</sub>",
                    )

        img = Draw.MolToImage(mol, size=size)
        img.save(filename)
