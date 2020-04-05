######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains the function to repair a PDB or a mmCIF file with pdbfixer
"""

import pdbfixer
import simtk.openmm.app

from HPC_Drug.auxiliary_functions import path

def _repair(input_file_name, file_type, output_file_name, add_H = False, ph = 7.0):
    """
    Private, it is called by the repair function

    repairs a PDB or mmCIF file with pdbfixer and returns the new file_name

    input_file_name :: string, the pdb or mmcif file to be repaired

    file_type :: string, can be cif or pdb

    output_file_name :: string, the name of the new structure file that will be created

    add_H :: bool, default False, if True pdbfixer will add hydrogens according to ph

    ph :: float, default 7.0, if add_H == True this is the pH value that will be used to add hydrogens
    """

    if file_type == 'cif':

        import HPC_Drug.PDB.biopython as bio

        #PART OF A PATCH TO ADRESS issue #195 on pdbfixer github
        TMP_dict = bio.mmcif2dict(file_name = input_file_name)

        patch_dict = {}

        for i in range(len(TMP_dict['_atom_site.label_asym_id'])):

            patch_dict[TMP_dict['_atom_site.label_asym_id'][i]] = TMP_dict['_atom_site.auth_asym_id'][i]
        #-----------------
        
        with open(input_file_name, 'r') as f:
            fixer = pdbfixer.PDBFixer(pdbxfile = f)
    
    elif file_type == 'pdb':

        with open(input_file_name, 'r') as f:
            fixer = pdbfixer.PDBFixer(pdbfile = f)

    else:
        raise TypeError(f"file_type must be cif or pdb not {file_type}")

    #calling pdbfixer methods
    fixer.findMissingResidues()

    fixer.findNonstandardResidues()

    fixer.replaceNonstandardResidues()

    #fixer.removeHeterogens(False)

    fixer.findMissingAtoms()

    fixer.addMissingAtoms()

    if add_H == True:

        fixer.addMissingHydrogens(ph)

    #fixer.addSolvent(fixer.topology.getUnitCellDimensions())

    
    if file_type == 'cif':

        with open(output_file_name, 'w') as f:
            simtk.openmm.app.pdbxfile.PDBxFile.writeFile(fixer.topology, fixer.positions, f, keepIds= True)


        #PART OF A PATCH TO ADRESS issue #195 on pdbfixer github
        TMP_dict = bio.mmcif2dict(file_name = output_file_name)

        for i in range(len(TMP_dict['_atom_site.label_asym_id'])):

            TMP_dict['_atom_site.auth_asym_id'][i] = patch_dict[TMP_dict['_atom_site.label_asym_id'][i]]
        
        bio.write_dict2mmcif(dictionary = TMP_dict, file_name = output_file_name)
        #---------------


    elif file_type == 'pdb':

        with open(output_file_name, 'w') as f:
            simtk.openmm.app.PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds = True)
         
    return output_file_name


def repair(Protein):
    """
    repairs a PDB or mmCIF file with pdbfixer and returns the new file_name

    This function calls _repair it is an interface to use a
    HPC_Drug.structures.protein.Protein instance on _repair in a simplified way

    Protein :: HPC_Drug.structures.protein.Protein instance

    return Protein
    """

    new_protein_file = _repair(input_file_name = Protein.pdb_file,
                            file_type = Protein.file_type,
                            output_file_name = f"{Protein.protein_id}_repaired.{Protein.file_type}",
                            add_H = False,
                            ph = 7.0)

    #get absolute path
    new_protein_file = path.absolute_filepath(path = new_protein_file)

    Protein.pdb_file = new_protein_file

    return Protein