"""
Copyright (c) 2020-2020 Maurice Karrenbrock

This software is open-source and is distributed under the
GNU Affero General Public License v3 (agpl v3) license

A copy of the license must be included with any copy of the program or part of it
"""

# This is the super class for Protein and Ligand

from HPC_Drug.PDB import biopython
from HPC_Drug.PDB import prody


class Structure(object):
    """Raises an Exception for not implemented structure type
    and contains some common methods"""

    def __init__(self):
                
        self.file_type = None

        self.pdb_file = None
        
        self.structure = None

        raise NotImplementedError('This structure was not implemented yet')

    def write(self, file_name = None, struct_type = 'biopython'):
        """
        This method writes self.structure on a {self.file_type} file (pdb, cif) using biopython (default) or prody (can only write pdb files)
        If no file_name is given self.pdb_file file will be overwritten otherwise a new file called file_name will be created and
        self.pdb_file will be updated with the new file_name

        file_name :: string, default self.pdb_file

        struct_type :: string, values: biopython (default), prody ; is the kind of structure in self.structure
        """

        if file_name == None:
            file_name = self.pdb_file

        if struct_type == 'prody':

            if self.file_type == 'pdb':

                prody.write_pdb(structure = self.structure, file_name = file_name)

            else:
                raise TypeError("Prody can only write pdb files not mmcif files, change Protein.file_type to pdb before you call this method")

        elif struct_type == 'biopython':

            biopython.write(structure = self.structure, file_type = self.file_type, file_name = file_name)

        self.pdb_file = file_name

    def update_structure(self, struct_type = "biopython"):
        """
        Parses the self.pdb_file file with the selected tool (biopython (default), mmcif2dict or prody) and updates self.structure
        prody can only parse pdb files, if you try to parse a cif with prody a TypeError will be raised
        mmcif2dict can only parse cif files, if you try to parse a pdb with mmcif2dict a TypeError will be raised

        structure_type :: string, values: biopython (default), prody, mmcif2dict
        """

        if struct_type == "biopython":

            self.structure = biopython.structure_factory(Protein = self)

        elif struct_type == "prody":

            if self.file_type != "pdb":
                raise TypeError("Prody can only parse pdb files")
            
            self.structure = prody.parse_pdb(file_name = self.pdb_file)

        elif struct_type == "mmcif2dict":

            self.structure = biopython.mmcif2dict(file_name = self.pdb_file)

        else:
            raise ValueError(f"Structure_type can be biopython, mmcif2dict or prody not {struct_type}")
