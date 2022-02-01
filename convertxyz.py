import os
from os import path
import sys
import getopt
import rdkit
from rdkit import Chem
import pandas as pd
import numpy as np
import subprocess
import argparse
from utils import b2bf
import time
import signal

def mk_xyz_from_smi(rn):
    command=f'obabel -ismi {rn}.smi -oxyz -O {rn}.xyz --gen3d'
    #pro=subprocess.Popen(command,shell=True,stdout=subprocess.PIPE,preexec_fn=os.setsid) #Using shell=True is a security hazard but I did it this way because I was being lazy
    pro=subprocess.Popen(command,shell=True) #Using shell=True is a security hazard but I did it this way because I was being lazy

    try:
        pro.wait(timeout=1)
    except subprocess.TimeoutExpired:
        pro.kill()


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
    

    '''The goal of this script is to be able to make a .xyz file for the specified smiles strings in a given spreadhseet
    '''
    #dir_path =os.path.dirname(os.path.realpath(__file__))
    orig_dir=os.path.dirname(os.path.abspath(inputfile))
    dir_path=os.path.dirname(os.path.abspath(inputfile))
    df=pd.read_excel(inputfile)
    #orig_dir=os.getcwd()

    #Get the indices
    index=list(df.index.values)
    react_type=['reactant','product']


    for i,rn in enumerate(index):

        for j,rt in enumerate(react_type):

            #Turn rn into a string with padded 0s at the beginning
            rn_str=str(rn+1)
            rn_0pad=rn_str.zfill(2) #in this case I have 2 digits total

            #Set the path where the smi file is
            if solv == 'solv':
                if solvationmethod != '':
                    rgpath=f'{dir_path}/{rn_0pad}/{rt}/{solv}/{functional}/{b2bf(basis)}/{solvationmethod}'
                else:
                    print('solvation method not specified. Examples incude -s CPCM')
                    return
            elif solv == 'gas':
                if solvationmethod== '':
                    rgpath=f'{dir_path}/{rn_0pad}/{rt}/{solv}/{functional}/{b2bf(basis)}'
                else:
                    print('Gas specified but an implicit solvation model was also specified. These are conflicting, please rectify this.')
                    return

            #Excecute the conversion from .smi file to .xyz file
            os.chdir(rgpath)
            mk_xyz_from_smi(rn_0pad)
            os.chdir(orig_dir)







if __name__=="__main__":
    tic=time.perf_counter()

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
    toc=time.perf_counter()

    print(f'Run time={toc-tic:0.4f}')
