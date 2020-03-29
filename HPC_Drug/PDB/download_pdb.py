"""
Copyright (c) 2020-2020 Maurice Karrenbrock

This software is open-source and is distributed under the
GNU Affero General Public License v3 (agpl v3) license

A copy of the license must be included with any copy of the program or part of it
"""

#This is the function used to download the pdb and mmcif files from the wwPDB

import Bio.PDB
import os


def download(protein_id, file_type = 'cif', pdir = None):
    """
    The function downloads a PDB or a mmCIF from wwwPDB in a selected directory
    the default directory is the working directory
    it returns the filename (str)
    
    protein_id :: string, it is the protein to download
    
    file_type :: string, it can be pdb or cif depending on the format required, default cif
    
    pdir :: string, default working directory, the directory where the file is saved
    
    return file_name , string
    
    raises a FileNotFoundError if the file is not downloaded correctly
    """
    

    if file_type == 'cif':
        file_type = 'mmCif'
    
    elif file_type != 'pdb' and file_type != 'cif' and file_type != 'mmCif':
        raise ValueError(f"Must be 'pdb' or 'cif' not {file_type}")
    
    if pdir == None:
        pdir = os.getcwd()

    pdbl = Bio.PDB.PDBList()
    file_name = pdbl.retrieve_pdb_file(protein_id, False, pdir, file_format = file_type, overwrite = True)
    
    if not os.path.exists(file_name):
        raise FileNotFoundError(f'Was not able to download the protein or to find {file_name}')
    
    else:
        return file_name