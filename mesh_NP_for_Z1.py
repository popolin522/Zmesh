# Authors: Popolin
# Purpose: Process LAMMPS data files and update bond information based on Blossom algorithm output.

import re
import os
import numpy as np

def read_lammps_data(file_path):
    """
    Reads a LAMMPS data file.

    Parameters:
        file_path (str): Path to the LAMMPS data file.

    Returns:
        tuple: Extracted header, masses, pair coefficients, bond coefficients, atoms, bonds, and velocities.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    header, masses, pair_coeffs, bond_coeffs, atoms, bonds, velocities = [], [], [], [], [], [], []
    current_section = None

    for line in lines:
        if "Atoms" in line:
            current_section = "atoms"
            continue
        elif "Bonds" in line:
            current_section = "bonds"
            continue
        elif "Masses" in line:
            current_section = "masses"
        elif "Pair Coeffs" in line:
            current_section = "pair_coeffs"
        elif "Bond Coeffs" in line:
            current_section = "bond_coeffs"
        elif "Velocities" in line:
            current_section = "velocities"

        if current_section == "atoms" and re.match(r'^\d+', line):
            atoms.append(line.strip().split())
        elif current_section == "bonds" and re.match(r'^\d+', line):
            bonds.append(line.strip().split())
        elif current_section == "masses":
            masses.append(line)
        elif current_section == "pair_coeffs":
            pair_coeffs.append(line)
        elif current_section == "bond_coeffs":
            bond_coeffs.append(line)
        elif current_section == "velocities" and re.match(r'^\d+', line):
            velocities.append(line.strip().split())
        else:
            header.append(line)

    if "Velocities\n" in header:
        header.remove("Velocities\n")

    return header, masses, pair_coeffs, bond_coeffs, atoms, bonds, velocities


def write_lammps_data(file_path, header, masses, pair_coeffs, bond_coeffs, atoms, bonds, velocities):
    """
    Writes data back into a LAMMPS format file.

    Parameters:
        file_path (str): Path to save the updated LAMMPS data file.
        header, masses, pair_coeffs, bond_coeffs, atoms, bonds, velocities: Data components to write.
    """
    with open(file_path, 'w') as file:
        file.writelines(header)
        file.writelines(masses)
        file.writelines(pair_coeffs)
        file.writelines(bond_coeffs)
        file.write('\nAtoms\n\n')
        for atom in atoms:
            file.write(' '.join(atom) + '\n')
        file.write('\nVelocities\n\n')
        for velocity in velocities:
            file.write(' '.join(velocity) + '\n')
        file.write('\nBonds\n\n')
        for bond in bonds:
            file.write(' '.join(bond) + '\n')


def add_new_bonds_based_on_blossom_algorithm(atoms, bonds, id_pair_lists, min_ref_id):
    """
    Adds new bonds based on the Blossom algorithm assignment for octahedron rigid bodies.

    Parameters:
        atoms (list): List of atoms data.
        bonds (list): List of existing bonds.
        id_pair_lists (ndarray): Pairs of atom IDs to bond.
        min_ref_id (int): Reference ID for tip atoms (user-provided value ).

    Returns:
        list: Updated bonds list.
    """
    atom_data = {}
    reference_tip_id_difference = {}

    for atom in atoms:
        atom_id = int(atom[0])
        mol_id = int(atom[1])
        atom_type = int(atom[2])
        if atom_type in [3, 5]:
            atom_data.setdefault(mol_id, []).append(atom_id)
            if atom_type == 5:
                if mol_id not in reference_tip_id_difference:
                    reference_tip_id_difference[mol_id] = atom_id
                else:
                    smallest_tip_id_of_octa = min(reference_tip_id_difference[mol_id], atom_id)
                    reference_tip_id_difference[mol_id] = smallest_tip_id_of_octa - min_ref_id

    bond_id = int(bonds[-1][0]) + 1 if bonds else 1
    for mol_id, atom_ids in atom_data.items():
        for id_tuple in id_pair_lists:
            atom_id_a = id_tuple[0] + reference_tip_id_difference[mol_id]
            atom_id_b = id_tuple[1] + reference_tip_id_difference[mol_id]
            bonds.append([str(bond_id), '1', str(atom_id_a), str(atom_id_b)])
            bond_id += 1

    return bonds


def process_lammps_file_with_mesh(input_file, output_file, id_pair_lists,min_ref_id):
    """
    Reads a LAMMPS file, adds new bonds based on Blossom algorithm, and writes the updated data to a new file.

    Parameters:
        input_file (str): Path to the input LAMMPS data file.
        output_file (str): Path to save the updated LAMMPS data file.
        id_pair_lists (ndarray): Pairs of atom IDs to bond.
    """
    header, masses, pair_coeffs, bond_coeffs, atoms, bonds, velocities = read_lammps_data(input_file)
    
    # Add new bonds
    bonds = add_new_bonds_based_on_blossom_algorithm(atoms, bonds, id_pair_lists,min_ref_id)
    
    # Update the number of bonds in the header
    header[4] = f"{len(bonds)} bonds\n"

    # Write the updated data
    write_lammps_data(output_file, header, masses, pair_coeffs, bond_coeffs, atoms, bonds, velocities)


if __name__ == "__main__":
    # Example usage
    save_dir = os.getcwd()
    id_pair_lists = np.load(os.path.join(save_dir, "meshed_rigid_body_pair_list.npy"))
    parent_path = ""
    input_LAMMPS_datafile = os.path.join(parent_path, 'simulation.data')
    output_LAMMPS_datafile = os.path.join(parent_path, 'simulation_with_meshed_body.data')
    min_ref_id = 38401 #user-provided value. Find it through ovitio (visual inspection). A typical choice would be the minimum id value in the reference rigid body. 
    process_lammps_file_with_mesh(input_LAMMPS_datafile, output_LAMMPS_datafile, id_pair_lists,min_ref_id)
