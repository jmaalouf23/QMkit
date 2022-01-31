import os
from os import path
import sys
import rdkit
import pandas as pd
import numpy as np


#The goal of this script is to be able to make a .smi file and .xyz file for the
#specified smiles scrints in a given spreadhseet
def mk_smi_file(smiles,dirnum):

    writepath = os.path.join(os.getcwd(),f'{dirnum}.smi')
    mode = 'a' if os.path.exists(writepath) else 'w'
    with open(writepath, mode) as f:
        f.truncate(0)
        f.write(smiles)

#Specify the excel file you want to read

#dir_path =os.path.dirname(os.path.realpath(__file__))
dir_path=os.getcwd()
df=pd.read_excel('PCET_RedoxPotentials.xlsx')
orig_dir=os.getcwd()


#Get the indices
index=list(df.index.values)
reagents=df['Reagent 1 (smiles)'].to_numpy()
products=df['Product 1 (smiles)'].to_numpy()

#Change this to specify the solvation in the directory
solv='gas' #set to either 'solv' or gas
functional='B3LYP'
basis='6-31G2dfp'
solvationmethod=None

for i,rn in enumerate(index):
    #pad index with leading 0s if necessary
    rn_str=str(rn+1)
    rn_0pad=rn_str.zfill(2) #in this case I have 5 digits total
    if solvationmethod != None:
        rgpath=f'{dir_path}/{rn_0pad}/reactant/{solv}/{functional}/{basis}/{solvationmethod}'
    else:
        rgpath=f'{dir_path}/{rn_0pad}/reactant/{solv}/{functional}/{basis}/'

    os.chdir(rgpath)
    mk_smi_file(reagents[i],rn_0pad)
    os.chdir(orig_dir)
    if solvationmethod != None:
        prpath=f'{dir_path}/{rn_0pad}/product/{solv}/{functional}/{basis}/{solvationmethod}'
    else:
        prpath=f'{dir_path}/{rn_0pad}/product/{solv}/{functional}/{basis}/'
    os.chdir(prpath)
    mk_smi_file(products[i],rn_0pad)
    os.chdir(orig_dir)
