import os
import csv

def read_atoms_section(lines):
    atoms = {}
    start_index = None
    end_index = None

    for i, line in enumerate(lines):
        if "[ atoms ]" in line:
            start_index = i + 1
        elif "[ bonds ]" in line:
            end_index = i
            break

    if start_index is not None and end_index is not None:
        for line in lines[start_index:end_index]:
            if line.strip() and not line.startswith(';'):
                data = line.split()
                number = int(data[0])
                atom_name = data[4]
                atoms[number] = atom_name

    return atoms

def print_assigned_atom_names(atoms):
    if not atoms:
        print("No assigned atom names found.")
    else:
        print("Assigned atom names for each number:")
        for number, atom_name in atoms.items():
            print(f"{number} = {atom_name}")

def assign_atom_names(itp_file):
    with open(itp_file, 'r') as file:
        lines = file.readlines()

    atoms = read_atoms_section(lines)
    print_assigned_atom_names(atoms)

def replace_numbers_with_atom_names(lines, atoms):
    atoms_section = False
    bonds_section = False
    angles_section = False
    dihedrals_section = False

    for i in range(len(lines)):
        line = lines[i]

        if line.strip().startswith('[ atoms ]'):
            atoms_section = True
            bonds_section = False
            angles_section = False
            dihedrals_section = False
        elif line.strip().startswith('[ bonds ]'):
            atoms_section = False
            bonds_section = True
            angles_section = False
            dihedrals_section = False
        elif line.strip().startswith('[ angles ]'):
            atoms_section = False
            bonds_section = False
            angles_section = True
            dihedrals_section = False
        elif line.strip().startswith('[ dihedrals ]'):
            atoms_section = False
            bonds_section = False
            angles_section = False
            dihedrals_section = True
        elif line.strip().startswith('['):
            atoms_section = False
            bonds_section = False
            angles_section = False
            dihedrals_section = False

        if atoms_section:
            if not line.startswith(';') and line.strip() and not line.strip().startswith('['):
                data = line.split()
                if len(data) > 4:
                    number = int(data[0])
                    if number in atoms:
                        data[4] = atoms[number]
                    lines[i] = ' '.join(data) + '\n'

        elif bonds_section or angles_section or dihedrals_section:
            if not line.startswith(';') and line.strip() and not line.strip().startswith('['):
                data = line.split()
                if bonds_section and len(data) >= 2:
                    for j in range(2):
                        if data[j].isdigit():
                            number = int(data[j])
                            if number in atoms:
                                data[j] = atoms[number]
                elif angles_section and len(data) >= 3:
                    for j in range(3):
                        if data[j].isdigit():
                            number = int(data[j])
                            if number in atoms:
                                data[j] = atoms[number]
                elif dihedrals_section and len(data) >= 4:
                    for j in range(4):
                        if data[j].isdigit():
                            number = int(data[j])
                            if number in atoms:
                                data[j] = atoms[number]

                lines[i] = ' '.join(data) + '\n'

    return lines

def save_changes(itp_file, lines):
    output_file = os.path.splitext(itp_file)[0] + '.txt'
    with open(output_file, 'w') as file:
        file.writelines(lines)

    print(f"New txt file has been created with atom names and saved to the output file: {output_file}")

def extract_sections(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    atoms_section = []
    bonds_section = []
    angles_section = []
    dihedrals_sections = []
    copy = False
    dihedrals_section = []

    for line in lines:
        if line.startswith('[ atoms ]'):
            copy = True
        elif line.startswith('[ bonds ]'):
            copy = False
            break
        elif copy and not line.startswith(';'):
            columns = line.strip().split()
            atoms_section.append('\t'.join([columns[4], columns[1], columns[6], columns[5]]))

    copy = False

    for line in lines:
        if line.startswith('[ bonds ]'):
            copy = True
        elif line.startswith('[ pairs ]') or line.startswith('[ angles ]'):
            copy = False
            break
        elif copy and not line.startswith(';'):
            columns = line.strip().split()
            bonds_section.append('\t'.join(columns[0:2] + columns[3:5]))

    copy = False

    for line in lines:
        if line.startswith('[ angles ]'):
            copy = True
        elif line.startswith('[ dihedrals ]'):
            copy = False
            break
        elif copy and not line.startswith(';'):
            columns = line.strip().split()
            angles_section.append('\t'.join(columns[0:3] + columns[4:6]))

    copy = False

    for line in lines:
        if line.startswith('[ dihedrals ]'):
            if dihedrals_section:
                dihedrals_sections.append(dihedrals_section)
                dihedrals_section = []
            copy = True
        elif line.startswith('[ exclusions ]'):
            copy = False
            break
        elif copy and not line.startswith(';'):
            columns = line.strip().split()
            dihedrals_section.append('\t'.join(columns[0:4] + columns[5:8]))

    if dihedrals_section:
        dihedrals_sections.append(dihedrals_section)

    new_filename = os.path.splitext(filename)[0] + '.rtp'

    with open(new_filename, 'w') as file:
        file.write('[ atoms ]\n')
        file.write('\n'.join(atoms_section))
        file.write('\n\n')

        file.write('[ bonds ]\n')
        file.write('\n'.join(bonds_section))
        file.write('\n\n')

        file.write('[ angles ]\n')
        file.write('\n'.join(angles_section))
        file.write('\n\n')

        for dihedrals_section in dihedrals_sections:
            file.write('[ dihedrals ]\n')
            file.write('\n'.join(dihedrals_section))
            file.write('\n\n')

    print(f'New RTP record for the selected file has been created: {new_filename}')

def list_rtp_files():
    rtp_files = [file for file in os.listdir(".") if file.endswith(".rtp")]
    return rtp_files

def print_file_list(file_list):
    print("Available RTP files:")
    for i, file in enumerate(file_list, 1):
        print(f"{i}. {file}")

def select_file(file_list):
    while True:
        try:
            choice = int(input("Please enter the number corresponding to the file you want to convert (or 0 to exit): "))
            if choice == 0:
                return None
            elif 1 <= choice <= len(file_list):
                return file_list[choice - 1]
            else:
                print("Invalid choice. Please select a valid number.")
        except ValueError:
            print("Invalid input. Please enter a number.")

def convert_rtp_to_csv(rtp_file):
    csv_file = os.path.splitext(rtp_file)[0] + ".csv"
    with open(rtp_file, 'r') as rtp_file, open(csv_file, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        for line in rtp_file:
            csv_writer.writerow(line.split())

def select_txt_file():
    txt_files = [f for f in os.listdir('.') if f.endswith('.txt')]

    if not txt_files:
        print('No .txt files found in the directory.')
        return

    print('Select the .txt file to create an RTP record from:')
    for i, file in enumerate(txt_files):
        print(f'{i+1}. {file}')

    try:
        option = int(input('Enter the number corresponding to the file: '))
        if option < 1 or option > len(txt_files):
            print('Invalid option.')
            return

        filename = txt_files[option - 1]
        extract_sections(filename)
    except ValueError:
        print('Invalid input.')

def main():
    itp_files = [file for file in os.listdir('.') if file.endswith('.itp')]

    if len(itp_files) == 0:
        print("No .itp files found in the current directory.")
        return

    print("Available .itp files:")
    for i, file in enumerate(itp_files):
        print(f"{i+1}. {file}")

    selection = input("Select the .itp file number to continue: ")

    if not selection.isdigit() or int(selection) <= 0 or int(selection) > len(itp_files):
        print("Invalid selection. Program terminated.")
        return

    selected_file = itp_files[int(selection) - 1]
    assign_atom_names(selected_file)
    proceed = input("Do you want to continue? (Yes/No): ")

    if proceed.lower() != 'yes':
        return

    atoms = read_atoms_section(open(selected_file, 'r').readlines())
    print_assigned_atom_names(atoms)

    while True:
        change = input("Would you like to make any changes? (Yes/No): ")

        if change.lower() != 'yes':
            break

        number = input("Type the number that you would like to change: ")
        if not number.isdigit() or int(number) not in atoms:
            print("Invalid number. Please select a correct number from the atom list.")
            continue

        atom_name = input("Type the atom name you would like to use for that number: ")
        atoms[int(number)] = atom_name

        print("Assigned atom names after change:")
        print_assigned_atom_names(atoms)

    lines = open(selected_file, 'r').readlines()
    modified_lines = replace_numbers_with_atom_names(lines, atoms)
    save_changes(selected_file, modified_lines)

    # Process the RTP file
    select_txt_file()
    
    # List and convert RTP files to CSV
    rtp_files = list_rtp_files()
    if rtp_files:
        print_file_list(rtp_files)
        chosen_file = select_file(rtp_files)
        if chosen_file:
            convert_rtp_to_csv(chosen_file)
            print(f"{chosen_file} successfully converted to CSV.")