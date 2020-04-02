"""
Copyright (c) 2020-2020 Maurice Karrenbrock

This software is open-source and is distributed under the
GNU Affero General Public License v3 (agpl v3) license

A copy of the license must be included with any copy of the program or part of it
"""

#Converts a mmcif in a pdb file

from HPC_Drug.PDB import biopython

def mmcif2pdb(Protein):
    """
    Takes a Protein instance and returns one
    if Protein.file_type == 'cif' converts the file in a pdb
    and updates Protein.file_type and Protein.pdb_file
    if Protein.file_type == 'pdb' does nothing.
    Otherwise raises TypeError
    
    Protein :: HPC_Drug.structures.protein.Protein instance
    
    return Protein
    """

    if Protein.file_type == 'pdb':
        pass

    elif Protein.file_type == 'cif':

        new_name = Protein.pdb_file.rsplit('.', 1)[0] + '.pdb'

        Protein.update_structure(struct_type = "biopython")

        biopython.write_mmcif(structure = Protein.structure, file_name = new_name)

        Protein.pdb_file = new_name

        Protein.file_type = 'pdb'

    else:
        raise TypeError(f"Protein.file_type must be 'pdb' or 'cif' not {Protein.file_type}")

    return Protein