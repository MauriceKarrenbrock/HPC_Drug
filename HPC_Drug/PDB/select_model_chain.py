######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

import HPC_Drug.PDB.biopython as _biopython

def select_model_chain(input_pdb, output_pdb, model=None, chain=None):
    """Select model and chain of a pdb or mmcif

    Parameters
    -------------
    input_pdb : str
    output_pdb : str
    model : int, default=None
        model to choose, if None none will be selected
    chain : str, default=None
        chain to choose, if None none will be selected

    Returns
    -----------
    output_pdb : str
    """

    input_file_type = str(input_pdb).strip().split('.')[-1]

    if input_file_type == 'pdb':
        structure = _biopython.parse_pdb(protein_id='aaaa', file_name=input_pdb)

    elif input_file_type == 'cif':
        structure = _biopython.parse_mmcif(protein_id='aaaa', file_name=input_pdb)

    else:
        raise ValueError(f'Only pdb and cif files are supported, not {input_file_type}')

    #select model
    if model is not None:
        structure = structure[model]

    #select chain
    if chain is not None:
        structure = structure[chain]

    _biopython.write(structure=structure,
                    file_type=str(output_pdb).strip().split('.')[-1],
                    file_name=output_pdb)

    return output_pdb
