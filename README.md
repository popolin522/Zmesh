# Zmesh
Python scripts to mesh rigid body (nanoparticles) for Z1+ algorithm. The scripts here can set up bonds between pairs of constituent atoms. The final output should be a processed data in LAMMPS data format ready to be used for Z1+algorithm

## How to use it?
1. Prepare a rigid body in LAMMPS data format (octa.data). This will be used as a reference
2. Install [NetworkX package](https://anaconda.org/anaconda/networkx#:~:text=To%20install%20this%20package%20run%20one%20of%20the) for set_mesh_pair.py
3. Run set_mesh_pair.py
4. Prepare the MD simulation file you with to mesh in LAMMPS data format (simulation.data), which is just a bunch of rigid bodies with a bunch of polymer.
5. Run mesh_NP_for_Z1.py

