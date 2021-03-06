######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains the class to run primadorac
"""

import os
import os.path
import shutil

from HPC_Drug.files_IO import read_file
from HPC_Drug.files_IO import write_on_files
from HPC_Drug.auxiliary_functions import path
from HPC_Drug.auxiliary_functions import run
from HPC_Drug.PDB import add_chain_id

class Primadorac(object):
    """
    It is a template class to run primadorac
    """

    def __init__(self, Protein, primadorac_path, ph = 7.0):
        """
        Protein :: Protein :: HPC_Drug.structures.protein.Protein instance
        with a valid Protein._ligands (Protein.get_ligand_list())

        primadorac_path :: string, the path to the primadorac executable (better if absolute)

        ph :: float, default 7.0, the ph at which the ligand will be protonated
        (at the moment primadorac does ONLY SUPPORT PH 7.0)
        """

        self.Protein = Protein

        #get the absolute path
        self.primadorac_path = path.absolute_programpath(program = primadorac_path)

        self.ph = ph


    def _run(self, string):
        """private"""
        
        run.subprocess_run(
            commands = string,
            shell = False,
            error_string = "primadorac failure"
        )

    def _rename_itp(self, file_to_search, ligand_resname = "LIG"):
        """
        private

        it is a patch because some old versions of primadorac do mess up the
        .itp file name
        """

        if not os.path.exists(file_to_search):

            new_name = f'{ligand_resname}.itp'

            if os.path.exists('file.itp'):

                os.rename('file.itp', new_name)
            
            elif os.path.exists('file.itp none'):
                
                os.rename('file.itp none', new_name)

            else:

                raise RuntimeError(f"Could not find any itp file created by primadorac\n\
                    Looked for {file_to_search}, file.itp and 'file.itp none'")

            return new_name

        else:

            return file_to_search

    def _edit_itp(self, ligand_resname, itp_file):
        """
        private

        primadorac itp call any lignd LIG i change it to the ligand_resname

        and removes the first 9 lines of the file (they make gromacs fail)
        """
        
        lines = read_file.read_file(file_name = itp_file)

        del lines[0:9]

        for i in range(len(lines)):

            lines[i] = lines[i].replace('LIG', ligand_resname)
            lines[i] = lines[i].replace('name-p', ligand_resname)
            lines[i] = lines[i].strip('\n')
            lines[i] = lines[i] + '\n'

        write_on_files.write_file(lines = lines, file_name = itp_file)

        return path.absolute_filepath(path = itp_file)

    
    def execute(self):
        """
        Run primadorac

        return Protein
        with updated .itp .prm .tpg files for any ligand in Protein._ligands
        """

        ligand_list = self.Protein.get_ligand_list()

        for i in range(len(ligand_list)):

            #primadorac can only work with files in it's workind directory, so if the file is somewhere else I copy it in the working directory
            try:

                current_pdb_file = shutil.copy(path.absolute_filepath(path = ligand_list[i].pdb_file), os.getcwd())

            except:

                current_pdb_file = os.path.relpath(ligand_list[i].pdb_file, os.getcwd())
            

            self._run(string = [self.primadorac_path, '-gp', current_pdb_file])



            prefix = ligand_list[i].pdb_file.strip().rsplit('.', 1)[0] + '-p'
            prefix = os.getcwd() + '/' + prefix.rsplit('/', 1)[-1]

            ligand_list[i].pdb_file = prefix + '.pdb'
            
            ligand_list[i].tpg_file = prefix + '.tpg'
            
            ligand_list[i].prm_file = prefix + '.prm'

            ligand_list[i].itp_file = self._rename_itp(
                file_to_search = prefix + '.itp',
                ligand_resname = ligand_list[i].resname
            )

            #any file_name is an absolute path


            ligand_list[i].itp_file = self._edit_itp(
                ligand_resname = ligand_list[i].resname,
                itp_file = ligand_list[i].itp_file
            )

            #adds the chain id in because primadorac removes it
            add_chain_id.add_chain_id(pdb_file = ligand_list[i].pdb_file, chain = self.Protein.chain)

        return self.Protein

            



