######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
The Gromacs input superclass
"""

import os

from HPC_Drug.files_IO import write_on_files
from HPC_Drug.files_IO import read_file
from HPC_Drug.auxiliary_functions import run
from HPC_Drug.MD.gromacs import gro2pdb

class GromacsInput(object):
    """Gromacs input template superclass
    
    The only public method is execute()"""

    def __init__(self,
                Protein = None,
                MD_program_path = 'gmx'):
        
        self.Protein = Protein

        #input filename is the cleaned gro file with both protein and ligand
        #but no solvent or trash HETATMS
        self.output_gro_file = os.getcwd() + f"{self.Protein.protein_id}_gromacs.gro"
        
        self.output_pdb_file = os.getcwd() +  f"{self.Protein.protein_id}_gromacs.pdb"
        
        self.MD_program_path = MD_program_path

        self.mdp_file = os.getcwd() + f"/{self.Protein.protein_id}_gromacs.mdp"
        
        self.command_string = []

        self.template = []

    

    def _write_template_on_file(self):
        """
        private
        
        Writes the object's template on self.mdp_file"""

        lines = ["\n".join(self.template)]

        write_on_files.write_file(lines = lines, file_name = self.output_pdb_file)
    

    def _interact_with_gromacs(self, string = None):
        """
        Private
        
        Interacts with gromacs running string
        with subprocess.run
        string must contain the gromacs path
        """

        print("Running Gromacs")

        run.subprocess_run(commands = string,
                        shell = False,
                        error_string = "Gromacs failure")


    def execute(self):

        self._write_template_on_file()

        self._interact_with_gromacs(string = self.command_string)

        try:
            #make a pdb file from the gro file
            self.Protein.pdb_file = gro2pdb.gro2pdb(gro_file = self.output_gro_file,
                                            pdb_file = self.output_pdb_file,
                                            chain = self.Protein.chain,
                                            gromacs_path = self.MD_program_path)

        except:
            pass

        self.Protein.gro_file = self.output_gro_file
 
        return self.Protein
