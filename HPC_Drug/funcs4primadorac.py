#contains classes and functions for primadorac

import structures
import pipeline_functions
import subprocess


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
                
                prefix = ligand_list[i].ligand_filename.rsplit('.', 1)[0]
                prefix = prefix + '-p'

                ligand_list[i].ligand_filename = prefix + '.pdb'
                ligand_list[i].topology_file = prefix + '.tpg'
                ligand_list[i].param_file = prefix + '.prm'

    return ligand_list