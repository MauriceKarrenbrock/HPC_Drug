######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

import os

from HPC_Drug.MD.gromacs import gromacs_input


class GromacsFirstOptimization(gromacs_input.GromacsInput):

    """Makes the first optimization of the given structure"""

    def __init__(self,
                Protein,
                MD_program_path = 'gmx'):

        super().__init__(Protein = Protein,
                        MD_program_path = MD_program_path)

        self.output_gro_file = os.getcwd() + f"{self.Protein.protein_id}_firstopt.gro"
        
        self.output_pdb_file = os.getcwd() +  f"{self.Protein.protein_id}_firstopt.pdb"

        self.mdp_file = os.getcwd() + f"/{self.Protein.protein_id}_firstopt.mdp"

        self.output_tpr_file = os.getcwd() +  f"{self.Protein.protein_id}_firstopt.tpr"

        #creates the .tpr and then optimizes the structure
        self.command_string = [

            f"{self.MD_program_path}", "grompp", "-f", f"{self.mdp_file}", "-c", f"{self.Protein.gro_file}", "-p", f"{self.Protein.top_file}", "-o", f"{self.output_tpr_file}", "-maxwarn", "100", "\n",
            
            f"{self.MD_program_path}", "mdrun", "-s", f"{self.output_tpr_file}", "-c", f"{self.output_gro_file}"
        ]

        self.template = ["integrator	= steep",
                        "nsteps		= 1000",
                        "emtol		= 100",
                        "emstep		= 0.01",
                        "nstxout 	= 1",
                        "nstenergy	= 1",
                        "rlist		= 1",
                        "coulombtype	= pme",
                        "vdw-type	= cut-off",
                        "rvdw		= 1.0",
                        "rcoulomb	= 1.0",
                        "constraints	= none"]
