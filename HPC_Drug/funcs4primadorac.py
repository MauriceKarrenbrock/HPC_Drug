#contains classes and functions for primadorac

from HPC_Drug import structures
from HPC_Drug import pipeline_functions

import subprocess
import os


def run_primadorac(ligand_list = None, primadorac_path = None, ph = 7.0):
    """Runs primadorac on any ligand given ad list
    or as simple Ligand istance
    
    ph choosing is not implemented yet
    will always be done for pH 7.0
    
    if a None or empty element is given as ligand_list None is returned"""

    if primadorac_path == None:
        raise Exception('Need a priamdorac path')

    if ligand_list == None or len(ligand_list) == 0:
        print("WARNING\nrun_primadorac received a None type list\
            \nwill keep going pretending nothing happened\nWARNING")
        return None

    print("Running primadorac")
    ligand_list = pipeline_functions.get_iterable(ligand_list)
    for i, ligand in enumerate(pipeline_functions.get_iterable(ligand_list)):
                #calls primadorac in order to get ligand's prm and tpg files
                #the -gp option means "don't optimize the structure but ad hydrogens for pH 7"
                #in the future primadorac should be able to add them at different pH too
                r = subprocess.run(f'bash {primadorac_path} -gp {ligand.ligand_filename}',
                                shell = True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                universal_newlines=True)

                print(r.stdout)
                print(r.stderr)

                if r.returncode != 0:
                    raise Exception(f"Primadorac failure\n{r.stdout}\n{r.stderr}")
                
                prefix = ligand_list[i].ligand_filename.rsplit('.', 1)[0]
                prefix = prefix + '-p'

                ligand_list[i].ligand_filename = prefix + '.pdb'
                ligand_list[i].topology_file = prefix + '.tpg'
                ligand_list[i].param_file = prefix + '.prm'
                ligand_list[i].itp_file = prefix + '.itp'

                #some old versions of primadorac do mess up the itp names
                if not os.path.exists(ligand_list[i].itp_file):
                    ligand_list[i].itp_file = rename_itp(ligand_resname = ligand_list[i].ligand_resname)

                #primadorac calls all ligands LIG inside the itp
                #this functions renames it to the right name
                ligand_list[i].itp_file = set_right_ligand_resname_in_itp(
                                                    ligand_resname = ligand_list[i].ligand_resname,
                                                    itp_file = ligand_list[i].itp_file)

    return ligand_list


def rename_itp(ligand_resname):
    """Renames the given itp file to ligand_resname.itp
    returns the new_name string"""

    new_name = f'{ligand_resname}.itp'

    if os.path.exists('file.itp'):

        os.rename('file.itp', new_name)
    
    elif os.path.exists('file.itp none'):
        
        os.rename('file.itp none', new_name)

    return new_name

def set_right_ligand_resname_in_itp(ligand_resname, itp_file):
    """primadorac itp call any lignd LIG i change it to the ligand_resname"""

    with open(itp_file, 'r') as f:
        text = f.readlines()

    with open(itp_file, 'w') as f:

        for line in text:

            line = line.replace('LIG', ligand_resname)
            f.write(f'{line}\n')

    return itp_file