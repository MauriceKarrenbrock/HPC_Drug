######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

from HPC_Drug.files_IO import write_on_files
from HPC_Drug.files_IO import read_file


def add_chain_id(pdb_file, chain = "A"):
    """
    This is a patch because orac and primadorac remove the chain id from
    pdb files and this confuses some pdb parsers (works on PDB files only)

    pdb_file :: string, the pdb file to edit

    chain :: string, default A, the chain id to add to the pdb_file

    returns nothing
    """

    chain = chain.upper().strip()

    lines = read_file.read_file(file_name = pdb_file)

    for i in range(len(lines)):

        if lines[i][0:4] == 'ATOM' or lines[i][0:6] == 'HETATM' or lines[i][0:3] == 'TER':

            lines[i] = lines[i][:20] + "{0:>2}".format(chain) + lines[i][22:]

            lines[i] = lines[i].strip('\n') + '\n'

    write_on_files.write_file(lines = lines, file_name = pdb_file)

