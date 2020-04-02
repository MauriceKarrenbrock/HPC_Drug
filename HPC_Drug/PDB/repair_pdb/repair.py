"""
Copyright (c) 2020-2020 Maurice Karrenbrock

This software is open-source and is distributed under the
GNU Affero General Public License v3 (agpl v3) license

A copy of the license must be included with any copy of the program or part of it
"""

#This file contains the function that repairs a pdb or mmcif file using the right repairing method

from HPC_Drug.PDB.repair_pdb import pdbfixer

def repair(Protein, repairing_method = "pdbfixer"):
    """
    This is a factory that returns a Protein with Protein.pdb_file updated to a repaired pdb or
    mmcif file using the right repiring method

    Protein :: HPC_Drug.structures.protein.Protein instance

    repairing_method :: string, default pdbfixer, it is the tool you want to repair the file with
    if you input a non existing tool will return NotImplementedError
    """

    if repairing_method == "pdbfixer":

        Protein = pdbfixer.repair(Protein = Protein)
        
    else:
        raise NotImplementedError

    return Protein