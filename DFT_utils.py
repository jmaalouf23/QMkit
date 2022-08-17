import os
import sys


def b2bf(basis_name):
    
    """
    Basis to basis folder: Converts a true basis function name to a name
    that is suitable for a folder in a directory. Removes *,(,)
    
    """

    if basis_name=='6-31G(2df,p)':
        return '6-31G2dfp'
    elif basis_name== '6-31++G**':
        return '6-31++Gstarstar'
    else:
        return basis_name

    
def bf2b(basis_name):
    
    """
    Basis folder to basis set name: Converts a folder name for a basis set into
    its true basis function name. 
    """

    if basis_name=='6-31G2dfp':
        return '6-31G(2df,p)'
    
    elif basis_name== '6-31++Gstarstar':
        return '6-31++G**'
    else:
        return basis_name