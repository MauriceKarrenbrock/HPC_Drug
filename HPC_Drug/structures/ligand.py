"""
Copyright (c) 2020-2020 Maurice Karrenbrock

This software is open-source and is distributed under the
GNU Affero General Public License v3 (agpl v3) license

A copy of the license must be included with any copy of the program or part of it
"""

#This is the Ligand class

from HPC_Drug.structures import structure


class Ligand(structure.Structure):
    """The Ligand class"""

    def __init__(self,
                resname = None,
                file_type = 'pdb',
                pdb_file = None,
                structure = None,
                resnum = None,
                itp_file = None,
                gro_file = None,
                top_file = None,
                tpg_file = None,
                prm_file = None):


        self.resname = resname
        
        self.file_type = file_type

        self.pdb_file = pdb_file
        if self.pdb_file == None:
            self.pdb_file = f"{self.resname}_lgand.{self.file_type}"

        self.structure = structure

        self.resnum = resnum
        if type(self.resnum) == str:
            self.resnum = int(self.resnum.strip())

        self.itp_file = itp_file

        self.gro_file = gro_file

        self.top_file = top_file

        self.tpg_file = tpg_file

        self.prm_file = prm_file
