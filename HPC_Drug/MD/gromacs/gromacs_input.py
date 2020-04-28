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
from HPC_Drug.auxiliary_functions import run
from HPC_Drug.auxiliary_functions import path
from HPC_Drug.MD.gromacs import gro2pdb

class GromacsInput(object):
    """Gromacs input template superclass
    
    The only public method is execute()"""

    def __init__(self,
                Protein = None,
                MD_program_path = 'gmx'):
        
        self.Protein = Protein
        
        self.output_gro_file = os.getcwd() + f"{self.Protein.protein_id}_gromacs.gro"
        
        self.output_pdb_file = os.getcwd() +  f"{self.Protein.protein_id}_gromacs.pdb"

        self.mdp_file = os.getcwd() + f"/{self.Protein.protein_id}_gromacs.mdp"
        
        self.MD_program_path = path.absolute_programpath(program = MD_program_path)
        
        self.command_string = []

        self.template = []

    
    def _gro2pdb(self, gro_file = None, pdb_file = None):
        """
        private
        """

        if gro_file is None:
            gro_file = self.output_gro_file

        if pdb_file is None:
            pdb_file = self.output_pdb_file

        pdb_file = gro2pdb.gro2pdb(gro_file = gro_file, pdb_file = pdb_file, chain = self.Protein.chain, gromacs_path = self.MD_program_path)

        return pdb_file

    def _pdb2gro(self, pdb_file = None, gro_file = None):
        """
        private
        """

        if gro_file is None:
            gro_file = self.output_gro_file

        if pdb_file is None:
            pdb_file = self.output_pdb_file

        gro_file = gro2pdb.pdb2gro(pdb_file = pdb_file, gro_file = gro_file, gromacs_path = self.MD_program_path)

        return gro_file

    def _write_template_on_file(self):
        """
        private
        
        Writes the object's template on self.mdp_file"""

        lines = ["\n".join(self.template)]

        write_on_files.write_file(lines = lines, file_name = self.mdp_file)
    

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
                        error_string = f"Gromacs failure\n{' '.join(self.command_string)}")


    def execute(self):

        self._write_template_on_file()

        self._interact_with_gromacs(string = self.command_string)

        #make a pdb file from the gro file
        self.Protein.pdb_file = self._gro2pdb()

        self.Protein.gro_file = self.output_gro_file
 
        return self.Protein
