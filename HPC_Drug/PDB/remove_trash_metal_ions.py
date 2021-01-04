######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains the function to remove trash metal ions
"""

from HPC_Drug import important_lists
from HPC_Drug.files_IO import read_file
from HPC_Drug.files_IO import write_on_files

def remove_trash_metal_ions(Protein, trash = important_lists.trash_ions):
    """This function removes unwanted metal ions
    that are still inside the structure after it went through prody
    selection (updates Protein.pdb_file)

    This is a brutal function I will need to do a better job
    
    Protein :: HPC_Drug.structures.protein.Protein instance

    Protein.file_type must be pdb or cif otherwise TypeError will be raised

    return Protein
    """

    def determine_pdb(line, trash = trash):

        #if it is a line containing atom coordinates
        if line[0:4] == 'ATOM' or line[0:6] == 'HETATM' or line[0:3] == 'TER':

            bool_output = line[17:20].strip().upper() not in trash

        #if it is not will not be kept anyway
        else:

            bool_output = False

        return bool_output

    def determine_cif(line, trash = trash):

        _line = line.split()

        #if it is a line containing atom coordinates
        if _line[0].strip() == 'ATOM' or _line[0].strip() == 'HETATM':

            bool_output = _line[5].strip().upper() not in trash

        #if it is not will be kept anyway
        else:

            bool_output = True

        return bool_output

    file_name = Protein.pdb_file

    lines = read_file.read_file(file_name = file_name)

    #keeping only the "good" lines
    if Protein.file_type == 'pdb':

        lines[:] = [x for x in lines if determine_pdb(x)]

        lines.append("end")

    elif Protein.file_type == 'cif':

        lines[:] = [x for x in lines if determine_cif(x)]

    else:
        raise TypeError(f"Protein.file_type must be pdb or cif not {Protein.file_type}") 

    for i in range(len(lines)):
        lines[i] = lines[i].strip() + '\n'

    write_on_files.write_file(lines = lines, file_name = file_name)

    Protein.pdb_file = file_name

    return Protein