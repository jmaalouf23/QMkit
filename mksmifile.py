import os
from os import path
import sys
import rdkit
import pandas as pd
import numpy as np
import argparse
from DFT_utils import b2bf

#The goal of this script is to be able to make a .smi file and .xyz file for the
#specified smiles scrints in a given spreadhseet



def mk_smi_file(smiles,dirnum):

    writepath = os.path.join(os.getcwd(),f'{dirnum}.smi')
    mode = 'a' if os.path.exists(writepath) else 'w'
    with open(writepath, mode) as f:
        f.truncate(0)
        f.write(smiles)



def main(argv):

    inputfile='' #Should be in the same directory as the numbered folders
    solv=''
    functional=''
    basis='' #Choosing this to be the default basis for the method. It is the one used by the QM9 dataset    
    solvationmethod=''

    inputfile=argv.i

    if argv.f is not None:
        functional=args.f
    if argv.b is not None:
        basis=args.b
    if argv.g is not None:
        solv=args.g
    if argv.s is not None:
        solvationmethod=argv.s
    if argv.info is not None:
        print('mksmifile.py -i <inputfile>  -g <gasorsolv> -f <functional>  -b <basis> -s <solvationmethod>')
    
    #Specify the excel file you want to read

    orig_dir=os.path.dirname(os.path.abspath(inputfile))
    dir_path=os.path.dirname(os.path.abspath(inputfile))
    df=pd.read_excel(inputfile)
    #orig_dir=os.getcwd()


    #Get the indices
    index=list(df.index.values)
    reagents=df['Reagent 1 (smiles)'].to_numpy()
    products=df['Product 1 (smiles)'].to_numpy()

    react_type=['reactant','product']
    react_lists=[reagents,products]
    
    for i,rn in enumerate(index):
        rgpath=''
        for j,rt in enumerate(react_type):
            #pad index with leading 0s if necessary
            rn_str=str(rn+1)
            rn_0pad=rn_str.zfill(2) #in this case I have 5 digits total
            if solv == 'solv':

                if solvationmethod != '':
                    rgpath=f'{dir_path}/{rn_0pad}/{rt}/{solv}/{functional}/{b2bf(basis)}/{solvationmethod}'
                else:
                    print('solvation method not specified. Examples incude -s CPCM')
                    return

            elif solv== 'gas':

                if solvationmethod=='':
                    rgpath=f'{dir_path}/{rn_0pad}/{rt}/{solv}/{functional}/{b2bf(basis)}/'
                else:
                    print('Gas specified but an implicit solvation model was also specified. These are conflicting, please rectify this.')
                    return
            os.chdir(rgpath)
            mk_smi_file(react_lists[j][i],rn_0pad)
            os.chdir(orig_dir)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True)
    parser.add_argument('-f', type=str)
    parser.add_argument('-b', type=str)
    parser.add_argument('-g', type=str)
    parser.add_argument('-s', type=str)
    parser.add_argument('--info', type=str)
    args = parser.parse_args()
    main(args)
    print("Done!")





