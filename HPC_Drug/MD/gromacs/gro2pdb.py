######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################


"""
gro2pdb and pdb2gro
"""

from HPC_Drug.PDB import add_chain_id
from HPC_Drug.auxiliary_functions import run


def gro2pdb(gro_file, pdb_file = None, chain = 'A', gromacs_path = 'gmx'):
    """converts a gro file in a pdb using gmx editconf and adds
    the chain identifier to the pdb (it is missing and it confuses many pdb parsers)

    gro_file :: string, the gro file to convert
    pdb_file :: string, optional, the name of the output file
    chain :: string, optional, default = A, it is the chain value that will be inserted in the pdb
    gromacs_path :: string, the gromacs executable default gmx

    all file names MUST be ABSOLUTE PATHS!!!

    returns a string of the pdb file name
    """

    if pdb_file == None:
        pdb_file = gro_file.rsplit('.', 1)[0].strip() + ".pdb"

    run.subprocess_run(commands = [f"{gromacs_path}", "editconf", "-f", f"{gro_file}", "-o", f"{pdb_file}"],
                    shell = False,
                    error_string = "Could not convert the gro file in pdb")

    add_chain_id.add_chain_id(pdb_file = pdb_file, chain = chain)

    return pdb_file



def pdb2gro(pdb_file, gro_file = None, gromacs_path = 'gmx'):
    """converts a pdb file in a gro using gmx editconf

    pdb_file :: string, the pdb file to convert
    gro_file :: string, optional, the name of the output file
    gromacs_path :: string, the gromacs executable default gmx

    all file names MUST be ABSOLUTE PATHS!!!

    returns a string of the gro file name
    """

    if gro_file == None:
        gro_file = pdb_file.rsplit('.', 1)[0].strip() + ".gro"

    run.subprocess_run(commands = [f"{gromacs_path}", "editconf", "-f", f"{pdb_file}", "-o", f"{gro_file}"],
                    shell = False,
                    error_string = "Could not convert the pdb file in gro")

    return gro_file


