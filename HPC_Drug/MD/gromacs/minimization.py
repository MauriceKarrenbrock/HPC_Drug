######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

import os

from HPC_Drug.MD.gromacs import gromacs_input


class GromacsMinimization(gromacs_input.GromacsInput):

    """Makes the first optimization of the given structure"""

    def __init__(self,
                Protein,
                MD_program_path = 'gmx'):

        super().__init__(Protein = Protein,
                        MD_program_path = MD_program_path)

        self.output_gro_file = os.getcwd() + "/" + f"{self.Protein.protein_id}_minimization.gro"
        
        self.output_pdb_file = os.getcwd() +  "/" + f"{self.Protein.protein_id}_minimization.pdb"

        self.mdp_file = os.getcwd() +  "/" + f"{self.Protein.protein_id}_minimization.mdp"

        self.output_tpr_file = os.getcwd() +  "/" + f"{self.Protein.protein_id}_minimization.tpr"

        #creates the .tpr and then optimizes the structure
        self.command_string = [

            [f"{self.MD_program_path}", "grompp", "-f", f"{self.mdp_file}", "-c", f"{self.Protein.pdb_file}", "-p", f"{self.Protein.top_file}", "-o", f"{self.output_tpr_file}", "-maxwarn", "100"],
            
            [f"{self.MD_program_path}", "mdrun", "-s", f"{self.output_tpr_file}", "-c", f"{self.output_gro_file}", '-ntmpi', '1']
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



class GromacsMinimizationOnlyLigand(gromacs_input.GromacsInput):

    """Makes the first optimization of all the ligands of a Protein"""

    def __init__(self,
                Protein,
                MD_program_path = 'gmx'):

        super().__init__(Protein = Protein,
                        MD_program_path = MD_program_path)

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

        self.mdp_file = os.getcwd() + f"/only_ligand_firstopt.mdp"

    def execute(self):

        if self.Protein.get_ligand_list() == []:
            return self.Protein

        self._write_template_on_file()

        Ligand = self.Protein.get_ligand_list()

        #creates the .tpr and then optimizes the structure
        #for any ligand
        for i in range(len(Ligand)):

            
            command_string = [
                [f"{self.MD_program_path}", "grompp", "-f", f"{self.mdp_file}", "-c", f"{Ligand[i].pdb_file}", "-p", f"{Ligand[i].top_file}", "-o", f"{os.getcwd()}/{Ligand[i].resname}_only_ligand.tpr", "-maxwarn", "100"],
                [f"{self.MD_program_path}", "mdrun", "-s", f"{os.getcwd()}/{Ligand[i].resname}_only_ligand.tpr", "-c", f"{os.getcwd()}/{Ligand[i].resname}_only_ligand.gro", '-ntmpi', '1']
            ]

            self._interact_with_gromacs(string = command_string)

            Ligand[i].gro_file = os.getcwd() + "/" + f"{Ligand[i].resname}_only_ligand.gro"

        return self.Protein


class MinimizeOnlyWaterBox(gromacs_input.GromacsInput):
    """
    Minimizes a box of solvent
    as deafault uses a copy of the one found in HPC_Drug.lib
    """

    def __init__(self,
                solvent_pdb,
                solvent_top,
                MD_program_path = "gmx"):

        self.solvent_pdb = solvent_pdb

        self.solvent_top = solvent_top

        self.MD_program_path = MD_program_path

        self.output_gro_file = os.getcwd() + "/" + f"only_water_{self.solvent_model}_minimization.gro"

        self.mdp_file = os.getcwd() + "/" + f"only_water_{self.solvent_model}_minimization.mdp"

        self.output_tpr_file = os.getcwd() + "/" + f"only_water_{self.solvent_model}_minimization.tpr"

        self.command_string = [
            [f"{self.MD_program_path}", "grompp", "-f", f"{self.mdp_file}", "-c", f"{self.solvent_pdb}", "-p", f"{self.solvent_top}", "-maxwarn", "100", "-o", f"{self.output_tpr_file}"],
            [f"{self.MD_program_path}", "mdrun", "-s", f"{self.output_tpr_file}", "-c", f"{self.output_gro_file}", '-ntmpi', '1']
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

    def execute(self):
        """
        Returns the updated gro and top files

        return gro_file, top_file
        """

        self._write_template_on_file()

        self._interact_with_gromacs()

        return self.output_gro_file, self.output_top_file
