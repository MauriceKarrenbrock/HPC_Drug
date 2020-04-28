######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################


"""
This file contains the classes and functions to get usable top and gro files starting from the given PDB and itp files
"""

import os
import subprocess

from HPC_Drug.MD.gromacs import gromacs_input
from HPC_Drug.files_IO import write_on_files
from HPC_Drug.files_IO import read_file
from HPC_Drug.PDB import merge_pdb

class GromacsMakeProteinGroTop(gromacs_input.GromacsInput):
    """
    Takes a protein instance and creates it's .top and .gro file
    for the Protein, ions and solvent BUT NOT for the organic ligands 
    THE ORGANIC LIGANDS SHOULD HAVE BEEN REMOVED BEFOREHAND FROM THE Protein.pdb_file
    """

    def __init__(self,
                Protein,
                solvent_model = '7',
                MD_program_path = "gmx"):

        super().__init__(self,
                Protein = Protein,
                MD_program_path = MD_program_path)

        self.solvent_model = solvent_model


        self.output_gro_file = os.getcwd() + f"/{self.Protein.protein_id}.gro"

        self.output_top_file = os.getcwd() + f"/{self.Protein.protein_id}.top"


        #protein_tpg_file is actually the protein model chosen through choices_file
        
        self.command_string = [f"{self.MD_program_path}", "pdb2gmx", "-f", f"{self.Protein.pdb_file}", "-o", f"{self.output_gro_file}", "-p", f"{self.output_top_file}"]

    
        self.template = [
            f"{self.Protein.tpg_file}",
            f"{self.solvent_model}"
                        ]
    
    def execute(self):

        """Takes a protein instance and returns one"""

        #write the file that will be given as choices for pdb2gmx (equivalent to pipeing the file with < on linux shell)
        write_on_files.write_file(lines = self.template, file_name = "choices.txt")

        #sometimes gromacs complains about some hydrogens, so in case I put the -ignh flag
        error_string = "Could not create the gromacs topology (.top) file"
        try:
            
            with open("choices.txt", "r") as std_in:
                
                r = subprocess.run(self.command_string,
                            shell = False,
                            stdin = std_in,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

            print(r.stdout)
            print(r.stderr)

            if r.returncode != 0:
                raise RuntimeError(error_string)

        except:

            ignh = self.command_string + ["-ignh"]

            with open("choices.txt", "r") as std_in:
                
                r = subprocess.run(ignh,
                            shell = False,
                            stdin = std_in,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

            print(r.stdout)
            print(r.stderr)

            if r.returncode != 0:
                raise RuntimeError(error_string)

        self.Protein.gro_file = self.output_gro_file
        self.Protein.top_file = self.output_top_file

        #update the pdb file in order to get the hydrogens
        self.Protein.pdb_file = self._gro2pdb(pdb_file = self.Protein.pdb_file)

        return self.Protein



class GromacsMakeJoinedProteinLigandTopGro(gromacs_input.GromacsInput):
    """
    Given a Protein instance containing some Ligand instances in Protein._ligands 
    and given that the protein already has a valid gro and top file
    adds the ligand(s) to the top file and to the Protein.gro_file file
    """

    def __init__(self,
                Protein,
                MD_program_path = "gmx"):

        super().__init__(self,
                Protein = Protein,
                MD_program_path = MD_program_path)

        self.output_gro_file = os.getcwd() + f"/{self.Protein.protein_id}_joined.gro"

        self.output_top_file = os.getcwd() + f"/{self.Protein.protein_id}_joined.top"

    def _edit_top_file(self):
        """
        Private
        
        Adds the needed #include and other informations
        to the protein top file in order
        to include the ligand
        """

        Ligand = self.Protein.get_ligand_list()

        top = read_file.read_file(file_name = self.Protein.top_file)

        itp_insertion_string = ''
        for lgand in Ligand:
            itp_insertion_string = itp_insertion_string + f'#include "{lgand.itp_file}"\n'
            compound_string = f'{lgand.resname}              1 \n'
            top.append(compound_string)

        for i in range(len(top)):
            if top[i].strip() == "" and top[i].strip()[0] == ";":
                pass
            elif top[i].strip().replace(" ", "").split(";")[0] == '[moleculetype]':
                top[i] = itp_insertion_string + '\n' + top[i]
                break

        write_on_files.write_file(lines = top, file_name = f'{self.Protein.protein_id}_joined.top')    

        self.Protein.top_file = f'{self.Protein.protein_id}_joined.top'
        
        return self.Protein


    def execute(self):

        #if there are no ligands do nothing
        if self.Protein.get_ligand_list() == []:
            return self.Protein
        
        #join the protein and ligand pdb
        self.Protein = merge_pdb.merge_pdb(Protein = self.Protein)

        #Create the joined gro file
        self.Protein.gro_file = self._pdb2gro(pdb_file = self.Protein.pdb_file)

        #edit and rename the top file
        self.Protein = self._edit_top_file()
        
        return self.Protein
