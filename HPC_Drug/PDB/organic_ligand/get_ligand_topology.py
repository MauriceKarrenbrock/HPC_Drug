######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains the function that uses the right tool to get
the topology of a given organic ligand (.itp .tpg .prm etc...)
in order to use it in a MD run
"""

from HPC_Drug.PDB.organic_ligand import primadorac
from HPC_Drug.auxiliary_functions import path

def get_topology(Protein, program_path, tool = "primadorac", ph = 7.0):
    """
    uses the right tool to get
    the topology of a given organic ligand (.itp .tpg .prm etc...)
    in order to use it in a MD run
    returns the updated protein

    Protein :: HPC_Drug.structures.protein.Protein instance
    with a valid Protein._ligands (Protein.get_ligand_list())
    value (if it is None or [] will return the protein untouched)

    program_path :: string, the absolute path to the tool's executable

    tool :: string, default primadorac, the tool to use to get the
    topology (.itp .tpg .prm etc...)

    ph :: float, default 7.0, the ph at which the ligand shall be added the missing hydrogens

    return Protein
    """

    if Protein.get_ligand_list() == None or Protein.get_ligand_list() == []:
        return Protein

    program_path = path.absolute_programpath(program = program_path)

    if tool == "primadorac":
        
        primadorac_obj = primadorac.Primadorac(Protein = Protein, primadorac_path = program_path, ph = ph)

        Protein = primadorac_obj.execute()

    else:

        raise NotImplementedError(f"{tool} is not an implemented tool")

    return Protein