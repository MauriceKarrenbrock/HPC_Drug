######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains the functions to parse and write files with biopython
"""


import Bio.PDB
import Bio.PDB.Entity
import Bio.PDB.MMCIF2Dict


######################################
# Parsers
#####################################

def parse_pdb(protein_id, file_name):
    """This function uses Biopython Bio.PDB.PDBParser to
    return a Biopython structure from a PDB file.

    protein_id :: string
    file_name :: string
    
    return structure """

    if type(protein_id) != str or type(file_name) != str:
        raise TypeError("Both protein_id and file_name must be strings")

    parser = Bio.PDB.PDBParser()

    structure = parser.get_structure(protein_id, file_name)

    return structure

def parse_mmcif(protein_id, file_name):
    """This function uses Biopython Bio.PDB.MMCIFParser to
    return a Biopython structure from a mmCIF file.

    protein_id :: string
    file_name :: string 
    
    return structure"""

    if type(protein_id) != str or type(file_name) != str:
        raise TypeError("Both protein_id and file_name must be strings")

    parser = Bio.PDB.MMCIFParser()

    structure = parser.get_structure(protein_id, file_name)

    return structure

def get_structure(file_name):
    """This is a function that uses the right
    parse_{} function depending on `file_name` suffix

    Parameters
    ------------
    file_name : str
        a pdb or mmcif file
    
    Returns
    ------------
    a biopython structure
    """
    
    file_suffix = str(file_name).split('.')[-1]

    if file_suffix == 'cif':
        return parse_mmcif(protein_id='aaaa', file_name=file_name)

    elif file_suffix in ('pdb', 'ent'):
        return parse_pdb(protein_id='aaaa', file_name=file_name)

    else:
        raise TypeError(f"file_name can only be 'cif' or 'pdb' not {file_suffix}")


def mmcif2dict(file_name):
    """Uses Bio.PDB.MMCIF2DICT to return a dictionary
    of the mmcif file"""

    if type(file_name) != str:
        raise TypeError("file_name must be a string")

    return Bio.PDB.MMCIF2Dict.MMCIF2Dict(file_name)


##################################################
# Writers
##################################################

def write_pdb(structure, file_name = "file.pdb"):
    """writes a pdb file when given a Biopython structure

    structure :: is instance Bio.PDB.Entity.Entity

    file_name :: string, default file.pdb

    returns nothing"""

    # Structure, Model, Chain, Residue
    if isinstance(structure, Bio.PDB.Entity.Entity):
        pass
    # List of Atoms
    elif hasattr(structure, '__iter__') and [x for x in structure if x.level == 'A']:
        pass
    else: # Some other weirdo object
        raise TypeError("Need a Bio.PDB.Entity instance like:\n"
                        "Structure, Model, Chain, Residue, list of Atoms.")

    if type(file_name) != str:
        raise TypeError("file_name must be a string")

    writer = Bio.PDB.PDBIO()

    writer.set_structure(structure)

    writer.save(file_name)


def write_mmcif(structure, file_name = "file.cif"):
    """writes a mmCIF file when given a Biopython structure

    structure :: is instance Bio.PDB.Entity.Entity

    file_name :: string, default file.cif

    returns nothing"""

    # Structure, Model, Chain, Residue
    if isinstance(structure, Bio.PDB.Entity.Entity):
        pass
    # List of Atoms
    elif hasattr(structure, '__iter__') and [x for x in structure if x.level == 'A']:
        pass
    else: # Some other weirdo object
        raise TypeError("Need a Bio.PDB.Entity instance like:\n"
                        "Structure, Model, Chain, Residue, list of Atoms.")

    if type(file_name) != str:
        raise TypeError("file_name must be a string")

    writer = Bio.PDB.MMCIFIO()

    writer.set_structure(structure)

    writer.save(file_name)

def write_dict2mmcif(dictionary, file_name = "file.cif"):
    """Writes a mmcif file starting from a dictionary
    obtained from mmcif2dict (that uses Bio.PDB.MMCIF2DICT )

    dictionary :: a dictionary containing all the mmcif infos, obtained with mmcif2dict

    returns nothing
    """

    if type(file_name) != str:
        raise TypeError("file_name must be a string")
    
    p = Bio.PDB.MMCIFIO()
    p.set_dict(dictionary)
    p.save(file_name)

def write(structure, file_type = "pdb", file_name = None):
    """This is a factory that writes the file
    given a structure or a mmcif2dict dictionary,
    the file type (pdb mmcif) and the output file name
    
    structure :: a Bio.PDB.Entity instance or a mmcif2dict dictionary
                 if you give a dictionary only file_type = 'cif' will be accepted

    file_type :: string, pdb or cif, default pdb

    file_name :: string, default file.{file_type}

    returns nothing
    """

    file_type = file_type.strip().lower()

    if file_name == None:
        file_name = f"file.{file_type}"

    if file_type == 'pdb':

        if type(structure) == dict:
            raise RuntimeError("If you give a dictionary as input file_type must be cif and not pdb")

        return write_pdb(structure = structure, file_name = file_name)

    elif file_type == 'cif':

        if type(structure) == dict:

            return write_dict2mmcif(dictionary = structure, file_name = file_name)

        else:

            return write_mmcif(structure = structure, file_name = file_name)

    else:

        raise TypeError(f"The file must be pdb or cif not {file_type}")