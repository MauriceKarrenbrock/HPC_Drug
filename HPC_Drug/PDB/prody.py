######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains the functions to write and parse pdb files with prody
"""

import prody

def parse_pdb(file_name):
    """Parses a PDB file with ProDy and returns a ProDy structure (prody.AtomGroup)

    file_name :: string

    returns structure
    """

    structure = prody.parsePDB(file_name)

    return structure

def write_pdb(structure, file_name = "file.pdb"):
    """
    Takes a Prody strucure prody.AtomGroup and writes it on a pdb file

    strucure :: prody.AtomGroup
    
    file_name :: string, default "file.pdb"

    returns nothing
    """

    prody.writePDB(file_name, structure)

def select(structure, string):
    """Uses the Prody select function with string as command

    structure :: prody.AtomGroup

    string :: string, this is the command that will be passed to prody select

    returns a new prody.AtomGroup
    """

    new_structure = structure.select(string)

    return new_structure

class ProdySelect(object):
    """This cass is a smart facade that implements some useful
    uses of HPC_Drug.PDB.prody.select
    """

    def __init__(self, structure):

        self._structure = structure

    def only_protein(self):
        """
        Returns a prody structure containing only the protein
        """

        return select(structure = self._structure, string = "protein")

    def protein_and_ions(self):
        """
        Returns a prody structure containing only the protein and the inorganic ions
        """

        try:

            return select(structure = self._structure, string = "protein or ion")

        except:

            return self.only_protein()

    def resname(self, resname):
        """
        Given a resname returns a Prody structure only containing any residue with
        that residue name

        resname :: string
        """

        resname = resname.strip().upper()

        return select(structure = self._structure, string = 'resname ' + resname)

    def resnum(self, resnum):
        """
        Given a resnum returns a Prody structure only containing the residue with
        that residue number

        resnum :: integer
        """

        return select(structure = self._structure, string = f"resnum {resnum}")