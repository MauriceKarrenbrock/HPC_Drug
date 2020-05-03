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

    it uses the force field inside Protein.tpg_file and must be the name of the .ff directory
    without the .ff suffix (case sensitive)

    solvent_model can e none spce tip3p etc (check the gromacs pdb2gmx -water flag documentation for extra info)
    """

    def __init__(self,
                Protein,
                solvent_model = "spce",
                MD_program_path = "gmx"):

        super().__init__(self,
                Protein = Protein,
                MD_program_path = MD_program_path)

        self.solvent_model = solvent_model


        self.output_gro_file = os.getcwd() + f"/{self.Protein.protein_id}.gro"

        self.output_top_file = os.getcwd() + f"/{self.Protein.protein_id}.top"


        #protein_tpg_file is actually the protein model chosen through choices_file
        
        self.command_string = [f"{self.MD_program_path}", "pdb2gmx",
            "-ff", f"{self.Protein.tpg_file}",
            "-water", f"{self.solvent_model}",
            "-f", f"{self.Protein.pdb_file}",
            "-o", f"{self.output_gro_file}",
            "-p", f"{self.output_top_file}"]

    
        self.template = []
    
    def execute(self):

        """Takes a protein instance and returns one"""

        #sometimes gromacs complains about some hydrogens, so in case I put the -ignh flag
        try:
            
            self._interact_with_gromacs(string = self.command_string)

        except:

            self._interact_with_gromacs(string = self.command_string + ["-ignh"])

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



class GromacsMakeOnlyLigandTopGro(gromacs_input.GromacsInput):

    """
    Makes the top and gro for any Ligand in a Protein instance
    """

    def __init__(self,
                Protein,
                MD_program_path = "gmx",
                solvent_model = "spce"):

        super.__init__(Protein = Protein,
                    MD_program_path = MD_program_path)

        self.solvent_model = solvent_model

    def _create_top_file(self):

        Ligand = self.Protein.get_ligand_list()

        for i in range(len(Ligand)):

            top_file = [
                '\n',
                ';protein ff (needed for the water itp file\n',
                f'#include  "{self.Protein.tpg_file}.ff/forcefield.itp" \n',
                '\n',
                ';ligand itp file'
                f'#include "{Ligand[i].itp_file}"\n',
                '\n',
                ';water itp',
                f'#include  "{self.Protein.tpg_file}.ff/{self.solvent_model}.itp" \n',
                '\n',
                "[ system ]",
                "; Name",
                "Protein in water",
                '\n',
                '[ molecules ]',
                '; Compound        #mols',
                f'{Ligand[i].resname}              1'
            ]
        
            write_on_files.write_file(lines = top_file, file_name = os.getcwd() + "/" + f"{Ligand[i].resname}_only_ligand.top")

            Ligand[i].top_file = os.getcwd() + f"{Ligand[i].resname}_only_ligand.top"
    
    def _create_gro_files(self):

        Ligand = self.Protein.get_ligand_list()

        for i in range(len(Ligand)):

            Ligand[i].gro_file = self._pdb2gro(pdb_file = Ligand[i].pdb_file,
                                            gro_file = os.getcwd() + "/" + f"{Ligand[i].resname}_only_ligand.gro")


    def execute(self):

        self._create_top_file()

        self._create_gro_files()

        return self.Protein

            
