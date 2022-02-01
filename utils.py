import os
import sys


def b2bf(basis_name):

    if basis_name=='6-31G(2df,p)':
        return '6-31G2dfp'
    elif basis_name== '6-31++G**':
        return '6-31++Gstarstar'
    else:
        return basis_name