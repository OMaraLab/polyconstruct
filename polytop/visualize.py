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
        self.atom_mapping = {}
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

        for atom in self.topology.atoms:
            atom_name = atom.atom_name
            # if the atom is a virtual atom (X), then set the element to H for rdKit 
            element = "H" if atom.atom_type == "X" else re.sub("[^a-zA-Z]", "", atom_name)
            new_atom = Chem.Atom(element)
            self.atom_mapping[atom.atom_id] = mol.AddAtom(new_atom)

        for bond in self.topology.bonds:
            try:
                if bond.order == 2:
                    bond_type = Chem.rdchem.BondType.DOUBLE
                elif bond.order == 3:
                    bond_type = Chem.rdchem.BondType.TRIPLE
                else:
                    bond_type = Chem.rdchem.BondType.SINGLE
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
            if not element == "C":
                if atom.is_virtual:
                    mol_atom.SetAtomicNum(0)
                    mol_atom.SetProp("RDKIT_ATOM_SYMBOL", "X")
                    mol_atom.SetFormalCharge(0)
                    atom_label = f"X<sub>{index}</sub>"
                else:
                    atom_label = element
                    count_h = sum([1 for neighbour in atom.bond_neighbours() if neighbour.element == "H"])
                    if count_h == 1:
                        atom_label += f"H"
                    elif count_h > 1:
                        atom_label += f"H<sub>{count_h}</sub>"
            mol_atom.SetProp("atomLabel", atom_label)
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
        virtual_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]
        img = Draw.MolToImage(mol, size=size, highlightAtoms=virtual_atoms )
        img.save(filename)
