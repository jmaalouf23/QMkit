# QMkit

• Generation of 3D coordinates and conformers from smiles strings using RDkit <br>
• Automated generation of quantum chemical calculation input and execution scripts.

### Key Functions

`get_xyz_from_smiles_string()` : Generates a text file with xyz coordinates of the molecule generated with the MMFF94 forcefield using RDkit.

`gaussian()` : Generates a gaussian submission script using an xyz text file (in the format generated from `get_xyz_from_smiles_string()`

`gen_gauss_sub_script()`: Can generate a slurm submission script to submit a gaussian calculation on a high performance computing environment.


