import os
from os import path
import sys
import getopt
import rdkit
from rdkit import Chem
import subprocess
import pandas as pd
import argparse
from utils import b2bf


def GetSpinMultiplicity(Mol, CheckMolProp = True):
    """Get spin multiplicity of a molecule. The spin multiplicity is either
    retrieved from 'SpinMultiplicity' molecule property or calculated from
    from the number of free radical electrons using Hund's rule of maximum
    multiplicity defined as 2S + 1 where S is the total electron spin. The
    total spin is 1/2 the number of free radical electrons in a molecule.

    Arguments:
    Mol (object): RDKit molecule object.
    CheckMolProp (bool): Check 'SpinMultiplicity' molecule property to
    retrieve spin multiplicity.

    Returns:
    int : Spin multiplicity.

    """

    Name = 'SpinMultiplicity'
    if (CheckMolProp and Mol.HasProp(Name)):
        return int(float(Mol.GetProp(Name)))

    # Calculate spin multiplicity using Hund's rule of maximum multiplicity...
    NumRadicalElectrons = 0
    for Atom in Mol.GetAtoms():
        NumRadicalElectrons += Atom.GetNumRadicalElectrons()

    TotalElectronicSpin = NumRadicalElectrons/2
    SpinMultiplicity = 2 * TotalElectronicSpin + 1

    return int(SpinMultiplicity)


def mkgauss_input_from_xyz(rn,smiles,filename,solvorgas='gas',solvmethod=None,solvent=None,functional='B3LYP',basis='6-31++G**'):
    writepath = os.path.join(os.getcwd(),f'{rn}.xyz')

    mode = 'r' if os.path.exists(writepath) else 'w'
    with open(writepath, mode) as f:
        xyz=f.readlines()

    writepath = os.path.join(os.getcwd(),f'{filename}.com')

    mode = 'a' if os.path.exists(writepath) else 'w'
    mol=Chem.MolFromSmiles(smiles)
    pc=Chem.rdmolops.GetFormalCharge(mol)
    #get multiplcity
    mult=GetSpinMultiplicity(mol)

    with open(writepath, mode) as f:
        f.truncate(0)
        f.write("%Mem=30GB\n")
        f.write(f"%chk={filename}.chk\n")
        f.write("%NProcShared=32\n")
        f.write(f"#n {functional}/{basis} Opt Freq\n")
        if solvorgas=='solv':
            f.write(f"SCRF=({solvmethod},solvent={solvent})\n\n")
            f.write('solventopt\n\n')
        else:
            f.write('\n\n')
        f.write(f'{pc} {mult}\n')
        for i,line in enumerate(xyz):
            if i>1:
                f.write(line)
        f.write('\n')
        if solvorgas=='solv':
            f.write('RADII=BONDI\n')
        f.write('\n')
        f.close()




# Actual executable script

def main(argv):
    inputfile=''
    outfile='' #For example like: 'reactant.com' but without the .com
    functional=''
    basis='' #CHoosing this to be the default basis for the method. It is the one used by the QM9 dataset
    solvorgas=''
    solvationmethod=''

    inputfile=argv.i
    outfile=argv.o
    
    if argv.f is not None:
        functional=args.f
    if argv.b is not None:
        basis=args.b
    if argv.g is not None:
        solvorgas=args.g
    if argv.s is not None:
        solvationmethod=argv.s
    if argv.info is not None:
        print('DFT_Script_gen.py -i <inputfile>  -o <outputfile> -f <functional> -g<gasorsolv> -b <basis> -s <solvationmethod>')
    

    orig_dir=os.path.dirname(os.path.abspath(inputfile))

    try:
        df=pd.read_csv(f'{inputfile}')
    except:
        df=pd.read_excel(f'{inputfile}')


    #Get reagents
#Get the indices

    index=list(df.index.values)
    reagents=df['Reagent 1 (smiles)'].to_numpy()
    products=df['Product 1 (smiles)'].to_numpy()
    try:
        solvent=df['solvent'].to_numpy()
    except:
        pass

    react_type=['reactant','product']
    react_lists=[reagents,products]

    for i,rn in enumerate(index):
        for j,rt in enumerate(react_type):
            #Turn rn into a string with padded 0s at the beginning
            rn_str=str(rn+1)
            rn_0pad=rn_str.zfill(2) #in this case I have 5 digits total
            if solvorgas=='solv':
                if solvationmethod !='':
                    reactantdir=os.path.join(orig_dir,f'{rn_0pad}',rt,solvorgas,functional,b2bf(basis),solvationmethod)
                    os.chdir(reactantdir)
                    mkgauss_input_from_xyz(rn_0pad,react_lists[j][i],f'{rt}',functional=functional, basis=basis,solvorgas=solvorgas,solvmethod=solvationmethod,solvent=solvent[i])
                else:
                    print('solvation method not specified. Examples incude -s CPCM')
                    return
   
            elif solvorgas=='gas':
                if solvationmethod=='':
                    reactantdir=os.path.join(orig_dir,f'{rn_0pad}',rt,solvorgas,functional,b2bf(basis))
                    os.chdir(reactantdir)
                    mkgauss_input_from_xyz(rn_0pad,react_lists[j][i],f'{rt}',functional=functional, basis=basis,solvorgas=solvorgas,solvmethod=solvationmethod,solvent=None)
                else:
                    print('Gas specified but an implicit solvation model was also specified. These are conflicting, please rectify this.')
                    return
            

    os.chdir(orig_dir)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True)
    parser.add_argument('-o', type=str)
    parser.add_argument('-f', type=str)
    parser.add_argument('-b', type=str)
    parser.add_argument('-g', type=str)
    parser.add_argument('-s', type=str)
    parser.add_argument('--info', type=str)
    args = parser.parse_args()
    main(args)
    print("Done!")
