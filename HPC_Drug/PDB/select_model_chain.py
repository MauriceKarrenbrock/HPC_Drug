######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

import Bio.PDB

def select_model_chain(Protein):
    """Takes a Protein instance containing the filename of a PDB or a mmcif
    Returns a Protein instance with an updated pdb or mmcif file
    using biopython
    selects only a chosen model and chain
    
    Protein.chain must be a string
    Protein.model must be an integer
    
    Protein :: HPC_Drug.structures.protein.Protein instance

    return Protein
    """    

    Protein.update_structure(struct_type = "biopython")

    #select model
    try:
    
        Protein.structure = Protein.structure[Protein.model]

    except KeyError:
        pass

    #select chain
    try:

        Protein.structure = Protein.structure[Protein.chain]

    except KeyError:
        pass

    #overwrites the protein file
    Protein.write(file_name = Protein.pdb_file, struct_type = "biopython")

    return Protein
