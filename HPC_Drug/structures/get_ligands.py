######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
The function in this file updates Protein._ligands
"""

from HPC_Drug.structures import ligand
from HPC_Drug.PDB import prody

def get_ligands(Protein, ligand_resnames_resnums):
    """
    Takes a HPC_Drug.structures.protein.Protein instance and a ligand_resnames_resnums
    and updates Protein._ligands with the newly created HPC_Drug.structures.ligand.Ligand
    instances

    Protein :: HPC_Drug.structures.protein.Protein, Protein.file_type must be pdb !!!

    ligand_resnames_resnums :: nested list of type [ ['ligand_resname', ligand_resnum], [...], ... ]
    if == None or == [] no HPC_Drug.structures.ligand.Ligand instance will be added to Protein._ligands

    return Protein
    """

    #removes old ligands if any
    Protein.clear_ligands()

    #There are no ligands to add
    if ligand_resnames_resnums == None or ligand_resnames_resnums == []:
        
        return Protein

    if Protein.file_type != 'pdb':
        raise TypeError(f"This function works only on PDB files, not {Protein.file_type}")

    #check if ligand_resnames_resnums is a nested list as it should be
    try:
        ligand_resnames_resnums[0][0]
    except:
        ligand_resnames_resnums = [ligand_resnames_resnums]


    Protein.update_structure(struct_type = "prody")

    prody_select = prody.ProdySelect(structure = Protein.structure)

    for i, ligand_name_number in enumerate(ligand_resnames_resnums):

        Ligand = ligand.Ligand(
            resname = ligand_name_number[0],
            pdb_file = None,
            structure = None,
            file_type = 'pdb',
            tpg_file = None,
            prm_file = None,
            resnum = ligand_name_number[1]
        )

        Ligand.structure = prody_select.resnum(resnum = Ligand.resnum)

        Ligand.write(file_name = f"{Ligand.resname}_lgand{i}.pdb", struct_type = "prody")

        Protein.add_ligand(Ligand = Ligand)

    return Protein

