######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains the functions to update existing ligands pdb files, resnums, and structures
"""

import copy

from HPC_Drug.structures import get_ligands
from HPC_Drug.PDB.structural_information import mmcif_header


def update_ligands(Protein, chain_model_selection = False):
    """
    If you have a HPC_Drug.structures.protein.Protein instance with some self._ligands
    but the files they referr to are obsolete, and maybe the resnums are obsolete too,
    if the Protein.pdb_file is a file containing (the protein and) the new ligand coordinates and resnums
    this is the function you want to use

    returns a HPC_Drug.structures.protein.Protein instance with updated self._ligands

    Protein :: HPC_Drug.structures.protein.Protein instance

    chain_model_selection :: bool, default False, if True checks only Protein.chain and Protein.model from the PDB file
    """

    ligand_resnames = []

    old_Protein = copy.deepcopy(Protein)

    for lig in Protein.get_ligand_list():

        ligand_resnames.append(lig.resname)

    Protein.update_structure(struct_type = "biopython")

    if chain_model_selection:
        protein_chain = Protein.chain
        protein_model = Protein.model

    else:
        protein_chain = protein_model = None

    Protein = get_ligands.get_ligands(Protein = Protein,
                        ligand_resnames_resnums = mmcif_header.get_ligand_resnum(structure = Protein.structure,
                                                                                ligand_resnames = ligand_resnames,
                                                                                protein_chain = protein_chain,
                                                                                protein_model = protein_model))

    #in order not to loose the info about itp, gro, etc... files
    for i in range(len(old_Protein.get_ligand_list())):

        old_Protein.get_ligand_list()[i].resnum = Protein.get_ligand_list()[i].resnum

        old_Protein.get_ligand_list()[i].pdb_file = Protein.get_ligand_list()[i].pdb_file

        if old_Protein.get_ligand_list()[i].resname.upper() != Protein.get_ligand_list()[i].resname.upper():

            raise ValueError("Something went wrong during the ligand update")

    return old_Protein





