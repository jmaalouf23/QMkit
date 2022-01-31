import os
from os import path
import sys
import getopt
import rdkit
from rdkit import Chem
import pandas as pd
import numpy as np
import subprocess

def mk_xyz_from_smi(rn):
    command=f'obabel -ismi {rn}.smi -oxyz -O {rn}.xyz --gen3d'
    subprocess.Popen(command,shell=True) #Using shell=True is a security hazard but I did it this way because I was being lazy

'''The goal of this script is to be able to make a .xyz file for the specified smiles strings in a given spreadhseet
'''
dir_path =os.path.dirname(os.path.realpath(__file__))
df=pd.read_excel('PCET_RedoxPotentials.xlsx')
orig_dir=os.getcwd()

#Get the indices
index=list(df.index.values)
reagents=df['Reagent 1 (smiles)'].to_numpy()
products=df['Product 1 (smiles)'].to_numpy()

#Change this to specify the solvation in the directory
solv='gas' # set this to either gas or solv
functional='B3LYP'
basis='6-31G2dfp'
solvationmethod=None
for i,rn in enumerate(index):

    #Turn rn into a string with padded 0s at the beginning
    rn_str=str(rn+1)
    rn_0pad=rn_str.zfill(2) #in this case I have 2 digits total

    #Set the path where the smi file is
    if solv == 'solv':
        rgpath=f'{dir_path}/{rn_0pad}/reactant/{solv}/{functional}/{basis}/{solvationmethod}'
    elif solv == 'gas':
        rgpath=f'{dir_path}/{rn_0pad}/reactant/{solv}/{functional}/{basis}'
    #Excecute the conversion from .smi file to .xyz file
    os.chdir(rgpath)
    mk_xyz_from_smi(rn_0pad)
    os.chdir(orig_dir)

    if solv == 'solv':
        prpath=f'{dir_path}/{rn_0pad}/product/{solv}/{functional}/{basis}/{solvationmethod}'
    elif solv == 'gas':
        prpath=f'{dir_path}/{rn_0pad}/product/{solv}/{functional}/{basis}'

    os.chdir(prpath)
    mk_xyz_from_smi(rn_0pad)
    os.chdir(orig_dir)
