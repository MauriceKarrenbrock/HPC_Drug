######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains a template to get the residues near a metallic ion, disulf bonds and organic ligand's renames and resnumbers
from a Protein instance and repair the PDB (or mmCIF) file
"""


from HPC_Drug.PDB.structural_information import mmcif_header
from HPC_Drug.PDB.structural_information import scan_structure
from HPC_Drug.PDB.repair_pdb import repair

class InfoRepair(object):
    """
    This is a template to get the residues near a metallic ion,
    disulf bonds and organic ligand's renames and resnumbers
    from a Protein instance
    and REPAIR the PDB (or mmCIF) file
    """

    def __init__(self, Protein, repairing_method = "pdbfixer"):

        self.Protein = Protein

        self.repairing_method = repairing_method

        self.organic_ligand_list = []

    def _parse_header(self):
        """private"""
        
        self.Protein, self.organic_ligand_list = mmcif_header.get_metalbinding_disulf_ligands(Protein = self.Protein)

    def _parse_structure(self):
        """private"""
        
        self.Protein, self.organic_ligand_list = scan_structure.get_metalbinding_disulf_ligands(Protein = self.Protein)


    def _repair(self):
        """private"""
        
        self.Protein = repair.repair(Protein = self.Protein,
                                    repairing_method = self.repairing_method)


    def _pdb(self):
        """private"""

        self._repair()

        self._parse_structure()
        

    def _cif(self):
        """private"""

        #tries to parse the header
        #if something goes wrong scans the structure
        try:

            self._parse_header()

            self._repair()

        except:

            self._pdb()


    def get_info_and_repair(self):
        """
        Returns Protein and organic_ligand_list
        
        Protein.pdb_file is a repaired PDB or mmCIF file

        return Protein, organic_ligand_list
        """

        if self.Protein.file_type == 'cif':

            self._cif()
        
        elif self.Protein.file_type == 'pdb':

            self._pdb()

        else:

            raise TypeError(f"Protein.file_type must be pdb or cif, not {self.Protein.file_type}")

        return self.Protein, self.organic_ligand_list