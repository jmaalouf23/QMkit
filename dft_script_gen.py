import os
from os import path
import rdkit
from rdkit import Chem
from dft_utils import b2bf, bf2b

ATOMIC_SYMBOLS = ['H', 'He','B','C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Br', 'I','Ar']


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



def gaussian(xyz_file,smiles,solvorgas='gas',solvmethod=None,solvent=None,functional='B3LYP',basis='6-31++G**',mem='180GB',mult=True,charge=True,pop_analysis=False,single_point=False,chk=False,n_proc_shared=48,file_extension='.com'):
    
    """ 
    Make a gaussian input script (usually a .com file) from an xyz file similar to the format produced by
    open babel. Takes a file in the output format of mk_xyz_from_smiles_string or mk_xyz_from_smi.

    Params:
    
    xyz_file       : Path to xyz file with coordinates. Do not include .xyz file extension
    solvorgas      : If 'gas' no solvation is included, if 'solv' then solvation will be included with method solvmethod
    solvent        : The solvent to be included in the calculation
    functional     : The functional to be used in the calculation if DFT is to be used
    basis          : the basis set to be used in the DFT or QM calculation
    mem            : The amount of memory to be used in the calculation
    mult           : Multiplicity to be included in the file. Default is to be calculated with GetSpinMultiplicity(). 
                     Other values can be specified for purposes such as population analysis and calulating cation or anion properties.
    charge         : The charge of the molecule in the calculation. Default is to be calculated with RDkit. 
                     Other values can be specified for purposes such as population analysis and calulating cation or anion properties.
    pop_analysis   : Whether or not to perform population analysis.
    single_point   : Whether or not to only perform single point calculation. If True, no geometry optimization is performed.
    chk            : Whether or not to include a check point file in the calculation
    n_proc_shared  : How many processors to be used by the gaussian calculation
    file_extension : The file extension used to save the file.

    Returns:
    
    gaussian input file containing functional, basis set, xyz coordinates, etc. of the same name as
    the xyz file

    """
    
    external_basis=['Def2SVPD','Def2TZVPD']
    
    
    writepath = os.path.join(os.getcwd(),f'{xyz_file}.xyz')
    short_filename=os.path.basename(os.path.normpath(xyz_file))
    
    mode = 'r' if os.path.exists(writepath) else 'w'
    with open(writepath, mode) as f:
        xyz=f.readlines()

    writepath = os.path.join(os.getcwd(),f'{xyz_file}{file_extension}')
    
    mode = 'a' if os.path.exists(writepath) else 'w'
    mol=Chem.MolFromSmiles(smiles)
    mol_withHs=Chem.AddHs(Chem.MolFromSmiles(smiles))
    #make list of atom lables
    atoms=[]
    for a in mol_withHs.GetAtoms():
        atoms.append(a.GetSymbol())
    atom_set=list(set(atoms))
    
    
    if charge and isinstance(charge,bool):
        pc=Chem.rdmolops.GetFormalCharge(mol)
    else:
        pc=charge
        
    #get multiplcity
    if mult and isinstance(mult,bool):
        mult=GetSpinMultiplicity(mol)
        

    with open(writepath, mode) as f:
        f.truncate(0)
        f.write(f"%Mem={mem}\n")
        if chk:
            f.write(f"%chk={short_filename}.chk\n")
            
        f.write(f"%NProcShared={n_proc_shared}\n")
        
        if single_point:
            sp=''
        else:
            sp='Opt=(maxcycles=120)' 
        if pop_analysis:
            p='Pop=Full iop(3/33=1)'
        else:
            p=''
            
        #Add logic for custom bais set
        if basis in external_basis:
            f.write(f"#p {functional}/Gen {p} {sp} Freq\n")
        else:
            f.write(f"#p {functional}/{bf2b(basis)} {p} {sp} Freq\n")
        
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

        
        #Add custom basis set info
        count=0
        if basis in external_basis:
            par_path=os.path.abspath(os.path.join(__file__, os.pardir))
            path= os.path.join(par_path,'basis_sets',f'{basis}.gbs')
            
            
            with open(path) as fp:
                Lines = fp.readlines()
                for line in Lines:
                    if count > 11:
                        curr_line_list=line.strip().split()
                        if len(curr_line_list)>0:
                            if curr_line_list[0] in ATOMIC_SYMBOLS and len(curr_line_list) <= 2:
                                curr_atom=curr_line_list[0]
                            if curr_atom in atom_set:
                                f.write(f'{line.strip()} \n')
                    count = count + 1
        else:
            f.write('\n')
            
        if solvorgas=='solv':
            f.write('\n')
            f.write('RADII=BONDI\n')
            f.write('\n')
        f.write('\n')  
        f.close()

        
        
        
def gen_gauss_sub_script(path,partition='xeon-p8', time='5-0:00:00',num_cpus=48,mem_per_cpu=4000,g16root='/home/gridsan/groups/manthiram_lab/gaussian',GAUSS_SCRDIR='/home/gridsan/groups/manthiram_lab/scratch/$SLURM_JOB_NAME-$SLURM_JOB_ID'):
    
    """
    Generates a slurm submission script for a gaussian calculation file (.com file). 
    Right now this is written for the MIT supercloud.
    
    Params
    
    path        : path where submission script will be saved. Include the entire path along with the filename.
                  Dont include the file extension, this is set to .sh by default. 
    partition   : #SBATCH --partition , will depend on system being used. 
    time        : #SBATCH --time , how  much time to allocate to the calculation before it is automatically canceled
    num_cpus    : #SBATCH -c , how many cps to be allocated to the task.
    mem_per_cpu : #SBATCH --mem-per-cpu , how much memory to use per CPU. Units of MBs.
    g16root     : Directory where gaussian is sotred on the system. Folder should be named g16, dont include this in the path.
    GAUSS_SCDIR : Scratch directory where temporary files will be written during the gaussian calculation.
    
    
    
    Returns
    
    .sh file that is able to execute a gaussian command for a given gaussian input script with the same name (usually a .com file)
    """
    
    short_filename=os.path.basename(os.path.normpath(path))
    writepath = os.path.join(os.getcwd(),f'{path}.sh')
    short_filename=os.path.basename(os.path.normpath(path))
    mode = 'a' if os.path.exists(writepath) else 'w'

    with open(writepath, mode) as f:
        f.truncate(0)
        f.write("#!/bin/bash -l\n")
        f.write(f"#SBATCH --partition={partition}\n")
        f.write(f"#SBATCH -J {short_filename}\n")
        f.write(f"#SBATCH -N 1\n")
        f.write(f"#SBATCH -c {num_cpus}\n")
        f.write(f"#SBATCH --time={time}\n")
        f.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
        f.write(f"#SBATCH --exclusive\n\n")
        
        f.write(f"export g16root={g16root}\n\n")
        f.write(f"export PATH=$g16root/g16/:$g16root/gv:$PATH\n\n")
        
        f.write(f"echo \"Gaussian PATH\"\n")
        f.write(f"which g16\n")
        f.write(f"input={short_filename}.com\n")
        
        f.write(f'echo "============================================================"\n')
        f.write(f'echo "Job ID : $SLURM_JOB_ID"\n')
        f.write(f'echo "Job Name : $SLURM_JOB_NAME"\n')
        f.write(f'echo "Starting on : $(date)"\n')
        f.write(f'echo "Running on node : $SLURMD_NODENAME"\n')
        f.write(f'echo "Current directory : $(pwd)"\n')

        f.write(f'echo "============================================================"\n\n')
        
        f.write(f'export GAUSS_SCRDIR={GAUSS_SCRDIR}\n\n')
        
        f.write(f'export GAUSS_SCRDIR\n')
        f.write(f'. $g16root/g16/bsd/g16.profile\n\n')
                
                
        f.write(f'echo "GAUSS_SCRDIR : $GAUSS_SCRDIR"\n')
        f.write(f'mkdir -p $GAUSS_SCRDIR\n')
        f.write(f'chmod 750 $GAUSS_SCRDIR\n')
        f.write(f'g16 $input\n\n')
        
        f.write(f'rm -rf $GAUSS_SCRDIR')
        
        f.write('\n')
        f.close()


