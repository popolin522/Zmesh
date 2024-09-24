# Zmesher
Python scripts to mesh rigid bodies (nanoparticles) for the Z1+ algorithm. These scripts create bonds between pairs of constituent atoms, and the final output is a processed LAMMPS data file ready for use in the Z1+ algorithm.

## How to use it?
1. Prepare a rigid body in LAMMPS data format (e.g., `octa.data`). This will serve as a reference.
2. Install the [NetworkX package](https://anaconda.org/anaconda/networkx#:~:text=To%20install%20this%20package%20run%20one%20of%20the) for use in `set_mesh_pair.py`.
3. Run `set_mesh_pair.py`.
4. Prepare the MD simulation data file that you wish to mesh in LAMMPS data format (e.g., `simulation.data`), which consists of multiple rigid bodies and polymers.
5. Identify a constituent atom in the reference data file from step (1) that will be used to establish connections in the MD simulation data. Document its ID and update the `min_ref_id` in the `mesh_NP_for_Z1.py` script accordingly.
6. Run `mesh_NP_for_Z1.py`.

## Technical details
### What are types 3 and 5 in the `mesh_NP_for_Z1.py` script?
In this setup, rigid bodies are made up of atoms of types 3 and 5. Each rigid body contains two constituent atoms of type 5, and the atom with the smaller ID is used as a reference. Users can modify this as needed.

### What is meant by "populating the connection"?
The scripts utilize the blossom algorithm (maximum matching) to pair constituent atoms, ensuring that (1) the distances between all pairs are minimized, and (2) all constituent atoms are paired. The next step is to replicate this pairing relationship across all rigid bodies in a large-scale simulation. To do this, a characteristic atom/bead is identified in a rigid body, and the "relative pairing relationship" established in `set_mesh_pair.py` is converted into an absolute pairing relationship (i.e., bonds in the LAMMPS data file).

### What is Z1+ algorithm?
[An awesome topological data analysis package for entangled polymer system](https://doi.org/10.1016/j.cpc.2022.108567).
