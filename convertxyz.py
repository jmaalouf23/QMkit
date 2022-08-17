import os
from os import path
import sys
import numpy as np
import subprocess
import time
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

"""
Author: Joseph Maalouf

Set of functions to 

• generate 3D coordinates from smiles strings and write text files with those coordinates
• extract 3D coordinates from the output files of QM or DFT calculations

"""

def mk_xyz_from_smiles_string(smiles: str ,filename: str,numConfs: int =1,conf_id=False,randomseed=0xf00d,return_mol=False,write_file=True):
    
    """
    Uses RDkit to create an xyz textfile from a smiles string. The output format of the xyz file is the same as what is made by open babel.
    The textfile is a .xyz format containing MMFF optimized conformers.
    Params
    
    smiles      : smiles string
    filename    : path, including name of file (without xyz extension), of the xyz file to be written.
    numConfs    : number of conformers to be generated by RDkit. Default is 1 but it is recommended to generate more and pick the lowest energy converged conformer.        
    conf_id     : of the generated conformers, conf_id species the geometry of which conformer to use to generate the xyz file. If not specified the conformer id of the lowest energy conformer is used.
    randomseed  : A random seed is used to generate conformers. For reproducability we set it to a default value.
    write_file  : If True,  a text file containing in the .xyz format containing MMFF optimized conformers with its name specified by the filename is written. 
    
    Returns
    if return_mol is True: returns the mol object of the lowest energy conformer
    """
    
    mol=Chem.MolFromSmiles(smiles)
    mol=Chem.AddHs(mol)
    
    try:
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=1,randomSeed=randomseed)
        res = AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=0)
    except:
        #If molecule cannot be emmbedded remove chiral stereochemistry and try again in case this was the issue.
        if '@' in smiles:
            smiles=smiles.replace('@','')
            print(f'CHIRAL STEROCHEMISTRY REMOVED, SMILES={smiles}')
            mol=Chem.MolFromSmiles(smiles)
            mol=Chem.AddHs(mol)
            cids = AllChem.EmbedMultipleConfs(mol, numConfs=numConfs,randomSeed=randomseed)
            res = AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=0)
            
        else:
            print('EMBEDDING FAILED FOR UNKNOWN REASON')
     
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=numConfs,randomSeed=randomseed)
    res = AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=0)
    
    if conf_id ==False: #If False the lowest energy converged conformer ID is used, otherwise the specified conf_id will be used
        
        res_converged=[v for i,v in enumerate(res) if res[i][0]==0]
        
        if len(res_converged) > 0:
            conf_id=res.index(min(res_converged, key = lambda t: t[1]))
        else:
            # If there are no converged conformers, try again but with useRandomCoords=True and using maxAttempts=400 instead of the default 200.
            # This will just increase the liklihood of getting a conveged conformer. This case will likely not have to be used.
            cids=AllChem.EmbedMultipleConfs(mol,numConfs=100,randomSeed=0xf00d,maxAttempts=400,useRandomCoords=True)
            res = AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=0)
            res_converged=[v for i,v in enumerate(res) if res[i][0]==0]
            if len(res_converged) > 0:
                conf_id=res.index(min(res_converged, key = lambda t: t[1]))
            else:
                conf_id=res.index(min(res, key = lambda t: t[1]))
    xyz=get_xyz_from_mol(mol,conf_id=conf_id) #gets an array of xyz values
    
    if write_file:
        write_xyz_from_xyz_arr(xyz,filename) #Writes xyz text file
        
    if return_mol:
        return ConfToMol(mol, conf_id)
        
    

def get_xyz_from_mol(mol,conf_id=0,randomseed=0xf00d):
    
    """
    Extracts a numpy array of coordinates from a molecules.
    Returns a `(N, 4)` numpy array of 3d coords of given rdkit molecule.
    N is the number of atoms in the molecule.
    Dim 0 : atomic symbol
    Dim 1 : x coordinates
    Dim 2 : y coordinates
    Dim 3 : z coordinates
    
    Params
    mol     : rdkit molecule to extract coordinates for
    conf_id : the conformer id number from which to extract 3D coordinates
    Returns
    -------
    Numpy ndarray of shape `(N, 3)` where `N = mol.GetNumAtoms()`.
    """

    xyz = np.zeros((mol.GetNumAtoms(), 4),dtype=object)
    
    if conf_id is not None:
        conf = mol.GetConformer(id=conf_id)
    else:
        conf = mol.GetConformer()
    atoms=mol.GetAtoms()
    for i in range(conf.GetNumAtoms()):
        position = conf.GetAtomPosition(i)
        atom=mol.GetAtomWithIdx(i)
        xyz[i,0] = atom.GetSymbol()
        xyz[i, 1] = position.x
        xyz[i, 2] = position.y
        xyz[i, 3] = position.z
    return (xyz)


def write_xyz_from_xyz_arr(xyz,filename,include_total_atoms=False):
    
    """Writes an XYZ text file from and xyz array given by get_xyz_from_mol.
    The ouput format exactly matches the format given by open babel.
    
    ----------
    xyz: numpy array of shape `(N , 4)` where `N = mol.GetNumAtoms(). The first
        column should be the atom label and the next 3 are the x,y,z coordinates
        
    filename: path to .xyz text file where information will be saved
    
    Returns
    -------
    .xyz text file containing atom lables and x,y,z coordinates.
    """    
    #might have to make filename jus the full path
    writepath = os.path.join(os.getcwd(),f'{filename}.xyz')
    mode = 'a' if os.path.exists(writepath) else 'w'
    with open(writepath, mode) as f:
        f.truncate(0)
        
        n,m=np.shape(xyz)
        if include_total_atoms:
            f.write(f'{n}\n')
        else:
            f.write('\n')
        f.write('\n')
        for i in range(n):
            f.write(f'{xyz[i,0]}{xyz[i,1]:>17,.6f}{xyz[i,2]:>15,.6f}{xyz[i,3]:>15,.6f}\n')
        f.close()
        
    
    
    
    
def combine_labels_and_xyz(labels : list, xyz : list) -> np.ndarray:
    
    """Combines a list of atom labels and list of lists of 3D coordinates
    into a single N x D+1 np array which can be used in the function
    write_xyz_from_xyz_arr()
    ----------
    labels: List of atom labels,  of size (N,). Each label should
    be a string
    
    xyz: List of lists of atom coordinates, each inner list contains the 3
    coordinates of atom i. ie [[0,0,0,],[1,1,1]] for 2 atoms with coordinates (0,0,0)
    and (1,1,1)
    
    Returns
    -------
    Returns a combined numpy array of size N x D+1 where the first column are the labels
    and the remaining columns are the cooridnates, much like a .xyz file
    """ 
    
    
    
    arr1,arr2= np.expand_dims(labels,1),np.array(xyz) #This unflattens the labels array 
    
    n,d=np.shape(np.array(arr2))
    arr=np.empty((n,d+1),dtype='object')
    arr[:,[0]]=arr1
    arr[:,1:]=arr2
    
    return arr  
    

def mk_xyz_from_smi(rn,software='obabel'):
    
    """
    Creates an XYZ file from a smiles string using open babel software. 
    
    """
    
    
    if software=='obabel':
        command=f'obabel -ismi {rn}.smi -oxyz -O {rn}.xyz --gen3d'
        pro=subprocess.Popen(command,shell=True) #Using shell=True is a security hazard. This should be updated in the future to avoid this.

        try:
            pro.wait(timeout=1)
        except subprocess.TimeoutExpired:
            pro.kill()
    elif software=='rdkit':
        
        #rdkit versions are now implemented in other functions, see mk_xyz_from_smiles_string()
        pass
    
def ConfToMol(mol, conf_id):
    conf = mol.GetConformer(conf_id)
    new_mol = Chem.Mol(mol)
    new_mol.RemoveAllConformers()
    new_mol.AddConformer(Chem.Conformer(conf))
    return new_mol

    