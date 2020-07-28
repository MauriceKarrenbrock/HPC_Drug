######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################


# This is the main file

from HPC_Drug import get_input
from HPC_Drug import pipelines
from HPC_Drug.files_IO import read_file, write_on_files
from HPC_Drug.PDB import merge_pdb

import sys
import os

def elaborate_autodoc_protein_file(file_name):
    
    pdb_lines = read_file.read_file(file_name = file_name)

    ligand_lines = []

    for i in range(len(pdb_lines)):


        if "REMARK" in pdb_lines[i]:
            break

        ligand_lines.append(pdb_lines[i])


    #remove the ligand part
    pdb_lines[:] = pdb_lines[i:]

    #atom type counter for ligand
    atom_count = {}


    #remove final useless stuff (like end lines) in pdb
    for i in range(len(pdb_lines) -1, -1, -1):

        if len(pdb_lines[i])  >= 6:
            if pdb_lines[i][0:6].strip() in ("ATOM", "HETATM"):

                last_residue_number = int(pdb_lines[i][22:26].strip())

                break

    pdb_lines[:] = pdb_lines[:i + 1]


    for i in range(len(ligand_lines)):

        tmp_line = ligand_lines[i].split()

        #rename the atom names counting them
        if tmp_line[2] in atom_count.keys():

            atom_count[tmp_line[2]] += 1

        else:

            atom_count[tmp_line[2]] = 1

        tmp_line[2] = f"{tmp_line[2]}{atom_count[tmp_line[2]]}"


        #ATOM      1  N   ARG A   3     -59.695  23.876  -8.802  1.00  0.00      AP1  N
        ligand_lines[i] = "{:<6}{:>5}{:>4}  {:>3}{:>2}{:>4d}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}\n".format(
                tmp_line[0],
                tmp_line[1],
                tmp_line[2],
                tmp_line[3],
                "A",
                last_residue_number + 1,
                float(tmp_line[5]),
                float(tmp_line[6]),
                float(tmp_line[7]),
                0.0,
                0.0
            )


    pdb_lines = pdb_lines + ligand_lines

    ligand_lines = None

    write_on_files.write_file(lines = pdb_lines, file_name = file_name[:-4] + "_new.pdb")

    pdb_lines = None



def main():
    try:
        input_file_name = sys.argv[1]
    except IndexError as err:
        raise IndexError("Needs a filename as input")

    try:
        input_variables = get_input.ParseInputFromFile(input_file_name).input_variables
    except FileNotFoundError as err:
        raise FileNotFoundError("Did not find the file")
    except ValueError as err:
        raise ValueError(err.args)



    elaborate_autodoc_protein_file(input_variables["filepath"])

    input_variables["filepath"] = input_variables["filepath"][0:-4] + "_new.pdb"


    pipeline = pipelines.choose_pipeline(**input_variables)

    pipeline.execute()

if __name__ == "__main__":

    main()