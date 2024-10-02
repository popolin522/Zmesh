# Authors: Po-An Lin @ Duke University @ Arya Lab
# Contact: poan.lin@duke.edu
# Purpose: Process LAMMPS data, compute distance matrices with PBC, and find minimized surface pairs using the Blossom algorithm.

import re
import os
import networkx as nx
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


def distance_matrix_pbc(xyz, box_lengths):
    """
    Computes the distance matrix with periodic boundary conditions (PBC).

    Parameters:
        xyz (numpy.ndarray): Array of shape (N, 3) containing the xyz coordinates of N atoms.
        box_lengths (numpy.ndarray): Array of shape (3,) containing the lengths of the simulation box in x, y, z directions.

    Returns:
        numpy.ndarray: A (N, N) distance matrix considering PBC.
    """
    N = xyz.shape[0]
    dist_matrix = np.zeros((N, N))
    
    for i in range(N):
        for j in range(i + 1, N):
            delta = xyz[i] - xyz[j]
            # Apply minimum image convention for periodic boundaries
            delta -= np.round(delta / box_lengths) * box_lengths
            dist = np.sqrt(np.sum(delta ** 2))
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist  # Symmetric matrix
    
    return dist_matrix


def find_minimized_surface_pair(coordinate, box_lengths=1000, save_dir=''):
    """
    Finds the minimized surface pair using the Blossom algorithm and saves the result.

    Parameters:
        coordinate (numpy.ndarray): Array of atom coordinates.
        box_lengths (numpy.ndarray or float): Lengths of the simulation box. Not an important thing here so I put in 1000 as a dummpy parameter. if the rigid body is crossing the periodic boundary, then be careful with box_lengths
        save_dir (str): Directory to save the output.

    Returns:
        np.ndarray: Pairs of atom IDs representing the minimized surface pair.
    """
    # Calculate distance matrix with periodic boundary conditions
    dist_matrix = distance_matrix_pbc(coordinate, box_lengths)

    # Create a graph and add edges with negative distances (to minimize total distance)
    G = nx.Graph()
    
    for i in range(len(coordinate)):
        for j in range(i + 1, len(coordinate)):
            dist = dist_matrix[i][j]
            G.add_edge(i, j, weight=-dist)  # Use negative distance to maximize matching

    # Find the minimum distance pairing using the Blossom algorithm
    matching = nx.max_weight_matching(G, maxcardinality=True)

    # Extract atom IDs
    id_array = atoms[:, 0].astype(int)

    # Convert index tuples from the matching to actual atom ID pairs
    id_pair_lists = np.array([tuple(id_array[idx] for idx in index_tuple) for index_tuple in matching])
    
    # Save the matching pairs to a file
    np.save(os.path.join(save_dir, "meshed_rigid_body_pair_list.npy"), id_pair_lists)
    
    return id_pair_lists


if __name__ == "__main__":
    # Example usage
    save_dir = os.getcwd()
    parent_path = "" 
    input_file = os.path.join(parent_path, 'octa.data')
    
    # Read atom data from LAMMPS file
    _, _, _, _, atoms, _, _ = read_lammps_data(input_file)
    
    # Extract coordinates (assumed to be in columns 3 to 5)
    coordinates = np.array(atoms)[:, 2:6].astype(float)
    
    # Find minimized surface pair and save result
    find_minimized_surface_pair(coordinates)
