import os
import re
import numpy as np

def list_rtp_files():
    """Lists all .rtp files in the current working directory."""
    files = [f for f in os.listdir() if f.endswith('.rtp')]
    return files

def select_rtp_file(files):
    """Prompts user to select an .rtp file."""
    print("Available .rtp files:")
    for i, file in enumerate(files):
        print(f"{i+1}. {file}")
    selection = int(input("Select the .rtp file (by number): ")) - 1
    return files[selection]

def read_rtp_file(file):
    """Reads the .rtp file and returns its content."""
    with open(file, 'r') as f:
        return f.readlines()

def parse_atoms_section(lines):
    """Parses the [ atoms ] section of the .rtp file."""
    atom_section = []
    inside_atoms = False
    for line in lines:
        if '[ atoms ]' in line:
            inside_atoms = True
        elif inside_atoms and line.strip() and not line.startswith('['):
            atom_section.append(line.strip().split())
        elif inside_atoms and (line.startswith('[') or not line.strip()):
            break
    return atom_section

def show_atoms(atom_section):
    """Displays available atoms to the user."""
    print("Available atoms:")
    for atom in atom_section:
        print(f"{atom[0]} (Charge: {atom[2]})")

def get_atoms_to_keep():
    """Asks user which atoms to keep."""
    atoms_to_keep = input("Enter the atom labels you want to keep (comma-separated): ")
    return [atom.strip() for atom in atoms_to_keep.split(',')]

def remove_unwanted_atoms(atom_section, atoms_to_keep):
    """Removes atoms that are not in the atoms_to_keep list."""
    return [atom for atom in atom_section if atom[0] in atoms_to_keep]

def neutralize_charges(atom_section, tolerance):
    """Attempts to neutralize the charges by adjusting within the given tolerance."""
    charges = np.array([float(atom[2]) for atom in atom_section])
    total_charge = np.sum(charges)
    
    if abs(total_charge) < 1e-6:
        return charges  # Already neutral
    
    adjustment_needed = -total_charge
    
    # Try adjusting each charge within the tolerance to achieve neutrality
    for i, charge in enumerate(charges):
        possible_adjustment = charge * tolerance
        if abs(possible_adjustment) >= abs(adjustment_needed):
            new_charges = charges.copy()
            new_charges[i] += adjustment_needed
            if abs(np.sum(new_charges)) < 1e-6:
                return new_charges

    return None  # Unable to neutralize within the given tolerance

def write_modified_rtp_file(file, lines, atom_section):
    """Writes the modified content back to the .rtp file."""
    with open(file, 'w') as f:
        inside_atoms = False
        atom_index = 0
        for line in lines:
            if '[ atoms ]' in line:
                inside_atoms = True
                f.write(line)
            elif inside_atoms and line.strip() and not line.startswith('['):
                if atom_index < len(atom_section):
                    f.write(f"{atom_section[atom_index][0]}\t{atom_section[atom_index][1]}\t{atom_section[atom_index][2]:.3f}\t{atom_section[atom_index][3]}\n")
                    atom_index += 1
            elif inside_atoms and (line.startswith('[') or not line.strip()):
                inside_atoms = False
                f.write(line)
            else:
                f.write(line)

def main():
    files = list_rtp_files()
    if not files:
        print("No .rtp files found in the current directory.")
        return

    rtp_file = select_rtp_file(files)
    lines = read_rtp_file(rtp_file)
    atom_section = parse_atoms_section(lines)

    show_atoms(atom_section)
    atoms_to_keep = get_atoms_to_keep()

    modified_atom_section = remove_unwanted_atoms(atom_section, atoms_to_keep)

    new_charges = neutralize_charges(modified_atom_section, 0.05)
    if new_charges is None:
        new_charges = neutralize_charges(modified_atom_section, 0.10)
        if new_charges is None:
            print("Unable to neutralize charges within 10% tolerance.")
            return

    for i, charge in enumerate(new_charges):
        modified_atom_section[i][2] = charge

    write_modified_rtp_file(rtp_file, lines, modified_atom_section)
    print(f"Modified .rtp file saved as {rtp_file}")