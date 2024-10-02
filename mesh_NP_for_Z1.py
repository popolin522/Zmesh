# Authors: Po-An Lin @ Duke University @ Arya Lab
# Contact: poan.lin@duke.edu
# Purpose: Process LAMMPS data files and update bond information based on Blossom algorithm output.

import re
import os
import numpy as np
import argparse
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


def add_new_bonds_based_on_blossom_algorithm(atoms, bonds, id_pair_lists, min_ref_id, 
                                             rigid_body_atom_types_list=[3,5], reference_atom_type=5, 
                                             do_selection_based_on_type=True):
    """
    Adds new bonds based on the Blossom algorithm assignment for octahedron rigid bodies.

    Parameters:
        atoms (list): List of atom data, where each atom entry is a list [atom_id, mol_id, atom_type].
        bonds (list): List of existing bonds.
        id_pair_lists (ndarray): Pairs of atom IDs to bond (output from Blossom algorithm).
        min_ref_id (int): Reference ID for tip atoms (user-provided value).
        rigid_body_atom_types_list (list): List of atom types representing rigid bodies.
        reference_atom_type (int): Atom type used for identifying reference atoms in rigid bodies.
        do_selection_based_on_type (bool): Whether to select based on atom type.

    Returns:
        bonds (list): Updated list of bonds with new bonds added.
        rigid_body_atom_id_list (list): List of atom IDs that belong to the rigid body.
    """
    atom_data = {}
    reference_tip_id_difference = {}
    rigid_body_atom_id_list = []
    
    for atom in atoms:
        atom_id = int(atom[0])
        mol_id = int(atom[1])
        atom_type = int(atom[2])

        # If the atom is part of the rigid body
        if atom_type in rigid_body_atom_types_list:
            rigid_body_atom_id_list.append(atom_id)
            atom_data.setdefault(mol_id, []).append(atom_id)
            
            if do_selection_based_on_type:
                # Handle selection based on the reference atom type
                if atom_type == reference_atom_type:
                    if mol_id not in reference_tip_id_difference:
                        reference_tip_id_difference[mol_id] = atom_id
                    else:
                        # Keep track of the smallest reference atom ID
                        smallest_tip_id_of_octa = min(reference_tip_id_difference[mol_id], atom_id)
                        reference_tip_id_difference[mol_id] = smallest_tip_id_of_octa - min_ref_id
            else:
                # If selection is not based on type, track the smallest atom_id for each rigid body
                if mol_id not in reference_tip_id_difference:
                    reference_tip_id_difference[mol_id] = atom_id
                else:
                    reference_tip_id_difference[mol_id] = min(reference_tip_id_difference[mol_id], atom_id)

    bond_id = int(bonds[-1][0]) + 1 if bonds else 1
    for mol_id, atom_ids in atom_data.items():
        for id_tuple in id_pair_lists:
            atom_id_a = id_tuple[0] + reference_tip_id_difference[mol_id]
            atom_id_b = id_tuple[1] + reference_tip_id_difference[mol_id]
            bonds.append([str(bond_id), '1', str(atom_id_a), str(atom_id_b)])
            bond_id += 1

    return bonds, rigid_body_atom_id_list

def delete_surface_polymer_bond(atoms, bonds, rigid_body_atom_id_list):
    """
    This function removes bonds between atoms in a rigid body and polymer atoms connecting to it.
    
    Parameters:
    - atoms: A list of atoms (not used in the current logic but retained for future expansion).
    - bonds: A list of bonds, where each bond contains the atom IDs.
    - rigid_body_atom_id_list: A list of atom IDs that belong to the rigid body.
    
    Returns:
    - Updated bonds list with bonds involving rigid body atoms and non-rigid body atoms removed.
    """
    
    # Filter out bonds where the first atom is in the rigid body and the second is not
    bonds = [bond for bond in bonds if not (int(bond[2]) in rigid_body_atom_id_list and int(bond[3]) not in rigid_body_atom_id_list)]
    
    return bonds


def process_lammps_file_with_mesh(input_file, output_file, id_pair_lists, min_ref_id, rigid_body_atom_types_list=[3,5], reference_atom_type=5, 
                                             do_selection_based_on_type=True, do_delete_surface_polymer_bond=True):
    """
    Reads a LAMMPS file, adds new bonds based on Blossom algorithm, and writes the updated data to a new file.

    Parameters:
        input_file (str): Path to the input LAMMPS data file.
        output_file (str): Path to save the updated LAMMPS data file.
        id_pair_lists (ndarray): Pairs of atom IDs to bond.
    """
    header, masses, pair_coeffs, bond_coeffs, atoms, bonds, velocities = read_lammps_data(input_file)
    
    # Add new bonds
    bonds, rigid_body_atom_id_list = add_new_bonds_based_on_blossom_algorithm(atoms, bonds, id_pair_lists, min_ref_id, 
                                             rigid_body_atom_types_list, reference_atom_type, 
                                             do_selection_based_on_type)
    
    #Delete the bond between surface bead and the polymer bead.
    if do_delete_surface_polymer_bond:
        bonds = delete_surface_polymer_bond(atoms, bonds, rigid_body_atom_id_list)

    # Update the number of bonds in the header
    header[4] = f"{len(bonds)} bonds\n"

    # Write the updated data
    write_lammps_data(output_file, header, masses, pair_coeffs, bond_coeffs, atoms, bonds, velocities)

def main():
    # Setup argument parser
    parser = argparse.ArgumentParser(description="Process LAMMPS data files and mesh rigid bodies with polymer grafts.")
    
    # Command-line arguments
    parser.add_argument('--save_dir', type=str, default=os.getcwd(), help='Directory to save/load files. Defaults to current working directory.')
    parser.add_argument('--parent_path', type=str, default='', help='Parent path for input/output data files.')
    parser.add_argument('--input_file', type=str, default='simulation_example.data', help='Input LAMMPS data file.')
    parser.add_argument('--output_file', type=str, default='simulation_with_meshed_body.data', help='Output LAMMPS data file.')
    parser.add_argument('--min_ref_id', type=int, required=True, help='Minimum reference atom ID for the rigid body.')
    
    # New arguments for rigid body atom types and reference type
    parser.add_argument('--rigid_body_atom_types', type=int, nargs='+', default=[3, 5], 
                        help='List of atom types representing rigid bodies. Defaults to [3, 5].')
    parser.add_argument('--reference_atom_type', type=int, default=5, 
                        help='Atom type used as the reference in rigid bodies. Defaults to 5.')
    parser.add_argument('--do_selection_based_on_type', type=bool, default=True, 
                        help='Flag to determine if the atom selection is based on type. Defaults to True. If False then it assume the reference id per rigid body is the one with the smallest atom id')
    
    # Argument for deleting surface polymer bonds
    parser.add_argument('--do_delete_surface_polymer_bond', type=bool, default=True, 
                        help='Flag to delete surface polymer bonds. Defaults to True.')

    # Parse arguments
    args = parser.parse_args()

    # Load the meshed rigid body pair list
    id_pair_lists = np.load(os.path.join(args.save_dir, "meshed_rigid_body_pair_list.npy"))
    
    # Full paths for input and output files
    input_LAMMPS_datafile = os.path.join(args.parent_path, args.input_file)
    output_LAMMPS_datafile = os.path.join(args.parent_path, args.output_file)

    # Call the function to process the LAMMPS file with the parsed arguments
    process_lammps_file_with_mesh(
        input_LAMMPS_datafile, 
        output_LAMMPS_datafile, 
        id_pair_lists, 
        args.min_ref_id, 
        args.rigid_body_atom_types, 
        args.reference_atom_type, 
        args.do_selection_based_on_type, 
        args.do_delete_surface_polymer_bond
    )

if __name__ == "__main__":
    main()

# if __name__ == "__main__":
#     # Example usage
#     save_dir = os.getcwd()
#     id_pair_lists = np.load(os.path.join(save_dir, "meshed_rigid_body_pair_list.npy"))
#     parent_path = ""
#     input_LAMMPS_datafile = os.path.join(parent_path, 'simulation.data')
#     output_LAMMPS_datafile = os.path.join(parent_path, 'simulation_with_meshed_body.data')
#     min_ref_id = 38401 #user-provided value. Find it through ovitio (visual inspection). A typical choice would be the minimum id value in the reference rigid body. 
#     process_lammps_file_with_mesh(input_LAMMPS_datafile, output_LAMMPS_datafile, id_pair_lists,min_ref_id)
