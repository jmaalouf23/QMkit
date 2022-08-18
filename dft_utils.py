import os
import sys
import cclib


"""
Useful functions used in gen_3D.py and dft_script_gen.py

"""

def b2bf(basis_name:str)->str:
    
    """
    Basis to basis folder: Converts a true basis function name to a name
    that is suitable for a folder in a directory. Removes *,(,)
    
    params:
    basis_name : name of a basis function that can be used directly in gaussian input script
    
    returns
    
    basis_name : updated basis name that can be used as a folder name in a directory
    
    """

    if basis_name=='6-31G(2df,p)':
        return '6-31G2dfp'
    elif basis_name== '6-31++G**':
        return '6-31++Gstarstar'
    else:
        return basis_name

    
def bf2b(basis_name:str)-> str:
    
    """
    Basis folder to basis set name: Converts a folder name for a basis set into
    its true basis function name. 
    
    params
    
    basis_name : basis name that can be used as a folder name in a directory
    
    returns
    
    basis_name : updated name of a basis function that can be used directly in gaussian input script
    
    """

    if basis_name=='6-31G2dfp':
        return '6-31G(2df,p)'
    
    elif basis_name== '6-31++Gstarstar':
        return '6-31++G**'
    else:
        return basis_name

def checkvibs(vibs,asint=False)-> bool:
    """ 
    checks the vibrational frequencies of a DFT calculation to make sure that they are all
    positive and thus there are no imaginary frequencies.

    Arguments:
    vibs (cclib vibrational frequencies): vibrational requencies obtained from gaussian output file using cclib

    Returns:
    True if all frequencies are positive. False otherwise. 

    """
    if asint:
        return int(all(vibs > 0))
    else:
        return all(vibs > 0)
    
    
def check(file,check_string, count=False):
    
    """ 
    Checks if a the specified string is contained within a given file. If count is not specified only presence is checked, otherwise if count is an integer
    then the function checks that the string is contained for the given number of counts.

    Arguments:
    file (file): any text file
    check_string (str): Check if file contains this string
    count: Number of times to check for check_string. If False then just checks to see if string is present.
    
    Returns:
    True if check string is present and count=False. If count=n, then returns true if string is contained for that specified number of times. 
    False if check_string is not present or if not present for the specified number of counts.

    """  
    
    with open(file) as f:
        if not count:
            return check_string in f.read()
        else: 
            return  f.read().count(check_string)==count




def ConfToMol(mol, conf_id):
    
    """
    Given an rdkit mol with set of conformers,
    returns the confomer with confomer id conf_id as
    a mol object (instead of conformer object).
    
    params
    
    mol :     rdkit mol object with conformers generated
    conf_id: id number of confer of mol to be selected
    
    returns
    
    new_mol: selected conformer of mol returned as an rdkit
             mol object 
    """
    
    conf = mol.GetConformer(conf_id)
    new_mol = Chem.Mol(mol)
    new_mol.RemoveAllConformers()
    new_mol.AddConformer(Chem.Conformer(conf))
    return new_mol



def mol_with_atom_index(mol, rev=False):
    
    """
    Given an rdkit mol, this function sets the atom map numer as the atom
    index so that when the molecule is displayed, the atom index is displayed in
    place of the atomic symbol for each atom. This is useful for debuggin.
    
    params:
    
    mol: rdkit mol
    rev: If false, replaces atomic symbols with atom index.
         If true, replaces atomic index with atomic symbols. 
         
    returns:
    
    mol: same instance of mol with updated atom map numbers.
    
    """
    if not rev:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())
        return mol
    else:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(atom.GetSymbol())
        return mol


def no_to_symbol(no:int)-> str:
    
    """
    Maps atomic number to atomic symbol
    
    params
    
    no: atomic number of molecule as integer
    
    returns
    
    atomic symbol corresponding the the atomic number given by 'no'
    
    """
    
    no_to_symbol_dict={
    1:'H',
    6:'C',
    7:'N',
    8:'O',
    9:'F',
    15:'P',
    16:'S',
    17:'Cl',
    35:'Br',
    53:'I'    
    }
    
    assert no in no_to_symbol_dict, "ATOMIC NUMBER NOT IN CURRENT DICTIONARY"
    
    
    return no_to_symbol_dict[no]