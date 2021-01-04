######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains the functions to merge PDB or mmCIF
They are useful wen you need to merge one or more organic ligands with a protein
"""

from HPC_Drug.files_IO import read_file
from HPC_Drug.files_IO import write_on_files
from HPC_Drug.auxiliary_functions import path


def merge_pdb(Protein):
    """
    Will put all the given ligands after the protein and update the ligand resnums
    this function is brutal and memory consuming I should do it better in the future

    both the protein and the ligands should be in PDB files (no check will be done)

    Protein :: HPC_Drug.structures.protein.Protein instance with a valid _ligands value

    return Protein with updated Protein.pdb_file
    """

    protein_file = read_file.read_file(file_name = Protein.pdb_file)

    ligands = Protein.get_ligand_list()

    #get the index of the line with the last ATOM HETATM or TER line
    #and get the resnum of this last residue
    for i in range(len(protein_file) - 1, -1 , -1):

        if protein_file[i][0:4] == 'ATOM' or protein_file[i][0:6] == 'HETATM' or protein_file[i][0:3] == 'TER':

            #some TER lines are non standard and don't contain the residue number
            try:
                residue_number =  int(protein_file[i][22:26].strip())
            except:
                residue_number =  int(protein_file[i-1][22:26].strip())

            index_protein_file = i + 1

            break


    #create the ligands list of strings
    ligand_file = []
    for j in range(len(ligands)):

        residue_number = residue_number + 1

        #update resnum
        ligands[j].resnum = residue_number

        tmp_ligand = read_file.read_file(file_name = ligands[j].pdb_file)

        #update the residue numbers in the file
        for k in range(len(tmp_ligand)):

            tmp_ligand[k] = tmp_ligand[k][:22] + "{0:>4}".format(residue_number) + tmp_ligand[k][26:]

        ligand_file = ligand_file + tmp_ligand


    #insert the ligands in the right place of the protein_file list
    protein_file[index_protein_file:index_protein_file] =  ligand_file

    #be sure to get the right formatting
    for i in range(len(protein_file)):

        protein_file[i] = protein_file[i].strip('\n') + '\n'

    #overwrite Protein.pdb_file
    write_on_files.write_file(lines = protein_file, file_name = f"{Protein.protein_id}_joined.pdb")

    Protein.pdb_file = path.absolute_filepath(path = f"{Protein.protein_id}_joined.pdb")

    return Protein

