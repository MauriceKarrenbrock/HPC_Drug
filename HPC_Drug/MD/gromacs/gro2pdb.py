######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################


"""
gro2pdb and pdb2gro
"""

from pathlib import Path
import mdtraj
from HPC_Drug.PDB import add_chain_id


def gro2pdb(gro_file, pdb_file = None, chain = 'A'):
    """converts a gro file in a pdb using mdtraj and adds
    the given chain identifier to the pdb

    gro_file :: string or path, the gro file to convert
    pdb_file :: string or path, optional, the name of the output file
    chain :: string, optional, default = A, it is the chain value that will be inserted in the pdb

    returns a string of the pdb file name
    """

    gro_file = Path(gro_file).resolve()

    if pdb_file == None:
        pdb_file = gro_file.with_suffix('.pdb')

    mdtraj.load(str(gro_file)).save(str(pdb_file), force_overwrite=True)

    add_chain_id.add_chain_id(pdb_file = pdb_file, chain = chain)

    return str(pdb_file)



def pdb2gro(pdb_file, gro_file = None):
    """converts a pdb file in a gro using mdtraj

    pdb_file :: string, the pdb file to convert
    gro_file :: string, optional, the name of the output file

    returns a string of the gro file name
    """

    pdb_file = Path(pdb_file).resolve()

    if gro_file == None:
        gro_file = pdb_file.with_suffix('.gro')

    mdtraj.load(str(pdb_file)).save(str(gro_file), force_overwrite=True)

    return str(gro_file)


