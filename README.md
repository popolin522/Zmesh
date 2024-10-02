# Zmesh
Python scripts to mesh rigid bodies (nanoparticles) with polymer grafts for the Z1+ algorithm. These scripts create bonds between pairs of surface constituent atoms of rigid body, and the final output is a processed LAMMPS data file ready for use in the Z1+ algorithm.
<img src="illustration.png" alt="" width="400"/>
## Why should we mesh the rigid body (nanoparticle) for Z1+ algorithm?
We want to trick the algorithm to account for the presence of nanoparticles during the entanglement analysis by "meshing" the surface of the rigid bodies. This involves creating bonds between pairs of atoms, where each atom connects to only one other atom, forming dumbbells. These dumbbells represent the confining surfaces of the nanoparticles. Since Z1+ keeps the terminal beads of these dumbbells fixed, the rigid obstacles (nanoparticles) remain stationary throughout the analysis.  This ensures that the polymer chains are prevented from crossing into the rigid obstacles during the minimization process.
## How to use it?
1. Prepare a rigid body in LAMMPS data format (e.g., `octa.data`) WITHOUT POLYMER. This will serve as a reference. If user only has data format with both polymer and rigid body, a LAMMPS script is provided to remove polymer (`delete.in.lmp`).
2. Install the [NetworkX package](https://anaconda.org/anaconda/networkx#:~:text=To%20install%20this%20package%20run%20one%20of%20the) for use in `set_mesh_pair.py`.
3. Run `set_mesh_pair.py`.
4. Prepare the MD simulation data file that you wish to mesh in LAMMPS data format (e.g., `simulation.data`), which consists of multiple rigid bodies and polymers.
5. Identify a constituent atom in the reference data file from step (1) that will be used to establish connections in the MD simulation data. Document its ID and update the `min_ref_id` in the `mesh_NP_for_Z1.py` script accordingly.
6. Run `mesh_NP_for_Z1.py`.

## Technical details
### What are types 3 and 5 in the `mesh_NP_for_Z1.py` script?
In this setup, rigid bodies are made up of atoms of types 3 and 5. Each rigid body contains two constituent atoms of type 5, and the atom with the smaller ID is used as a reference. Users can modify this as needed.

### What is meant by "populating the connection"?
The scripts utilize the blossom algorithm (maximum matching) to pair constituent atoms, ensuring that (1) the distances between all pairs are minimized, and (2) all constituent atoms are paired. The next step is to broadcast this pairing relationship across all rigid bodies in a large-scale simulation. To do this, user should identify a characteristic atom/bead in a rigid body in `octa.data`. Then, the "relative pairing relationship" established in `set_mesh_pair.py` is converted into an absolute pairing relationship (i.e., bonds in the LAMMPS data file) by refering to that characteristic atom/bead.

### What is Z1+ algorithm?
[An awesome topological data analysis package for entangled polymer system](https://doi.org/10.1016/j.cpc.2022.108567).

### Why cleaving bonds between constituent atom and polymer?
For a polymer-grafted system, Zmesh will automatically remove the bond between the surface constituent atom and the polymer to prevent generating spurious statistics in Z1+. For instance, if we don't remove such a bond, Z1+ may mistakenly interpret a chain with beads N=5 grafted onto paired constituent atoms as N=7.
### what is polymer-grafted system and why do we care?
Polymer-grafted nanoparticle is an interesting system to study because they combine the mechanical properties / high processeibility of polymers and unique functionalities provided by nanoparticles.
### other assumptions made in this code
I assume each rigid body contains an even number of atoms, and that the atom ID pattern is consistent across all rigid bodies (i.e., the IDs of surface atoms were not randomized for different nanoparticles. I mean, who does that?)
