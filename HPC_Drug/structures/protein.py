######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains the Protein class
"""

from HPC_Drug.structures import structure

class Protein(structure.Structure):
    """The Protein class"""

    def __init__(self,
                protein_id = None,
                pdb_file = None,
                structure = None,
                substitutions_dict = None,
                sulf_bonds = None,
                seqres = None,
                file_type = 'cif',
                model = 0,
                chain = 'A',
                cys_dict = None,
                gro_file = None,
                top_file = None,
                tpg_file = None,
                prm_file = None):
                
        if protein_id == None:
            raise ValueError('need a protein id')
        else:
            self.protein_id = protein_id.strip()

        self.model = model
        if type(self.model) == str:
            self.model = int(self.model)

        self.chain = chain.strip().upper()

        self.file_type = file_type.strip().lower()

        self.pdb_file = pdb_file

        
        if self.pdb_file == None:
            self.pdb_file = self.protein_id + '.' + self.file_type
        
        self.structure = structure

        self.substitutions_dict = substitutions_dict
        self.sulf_bonds = sulf_bonds
        self.seqres = seqres

        self.cys_dict = cys_dict

        self.gro_file = gro_file
        self.top_file = top_file

        self.tpg_file = tpg_file
        self.prm_file = prm_file

        self._ligands = []

    def add_ligand(self, Ligand):
        """adds a Ligand instance to self._ligands"""

        self._ligands.append(Ligand)

    def clear_ligands(self):
        """Clears ALL the ligands stored in self._ligands"""

        self._ligands = []

    def update_ligands(self, ligands):
        """
        Takes an iterable (list, tuple, etc...) containing HPC_Drug.structures.ligand.Ligand
        instances and ovewrites self._ligand with this new ones (any information about the old ones will be lost)
        
        ligands :: iterable containing the new HPC_Drug.structures.ligand.Ligand instances
        """

        self.clear_ligands()

        for lig in ligands:
            self.add_ligand(Ligand = lig)
            

    def get_ligand_list(self):
        """returns the list of ligands (self._ligands) already stored
        it is a pointer to it, not a copy, so pay attention"""

        return self._ligands