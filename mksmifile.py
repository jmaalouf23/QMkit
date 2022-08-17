import os
from os import path

#consider renaming make_smiles_file
def mk_smi_file(smiles: str , save_path: str): 
    
    """
    Given a smile, makes a file containing the smiles with extension .smi
    
    Params
    smiles    : smiles string to be written in file
    save_path : directory where file will be written. Do not include file extension
    
    """
    
    writepath = os.path.join(os.getcwd(),f'{save_path}.smi')
    mode = 'a' if os.path.exists(writepath) else 'w'
    
    with open(writepath, mode) as f:
        f.truncate(0)
        f.write(smiles)





