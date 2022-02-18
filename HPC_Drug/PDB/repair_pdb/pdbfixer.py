######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
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

def repair(input_file_name, output_file_name, add_H=False, ph=7.0):
    """
    Repairs a pdb file or a mmcif/pdbx

    repairs a PDB or mmCIF file with pdbfixer and returns the new file name
    will add missing atoms and residues, replace non standard atom and residue names
    with standard ones and if `add_H`=True adds hydrogens

    Parameters
    ------------
    input_file_name : str
        the pdb or mmcif file to be repaired
    output_file_name : str
        the repaired pdb or mmcif
    add_H : bool, default=False
        if True pdbfixer will add hydrogens according to ph
    ph : float, default=7.0
        if add_H == True this is the pH value that will be used to add hydrogens

    Returns
    ---------
    output_file_name
    """

    input_file_type = str(input_file_name).strip().split('.')[-1]

    output_file_type = str(output_file_name).strip().split('.')[-1]

    if input_file_type == 'cif':
        
        with open(input_file_name, 'r') as f:
            fixer = pdbfixer.PDBFixer(pdbxfile = f)
    
    elif input_file_type == 'pdb':

        with open(input_file_name, 'r') as f:
            fixer = pdbfixer.PDBFixer(pdbfile = f)

    else:
        raise TypeError(f"file type must be cif or pdb not {input_file_type}")

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

    
    if output_file_type == 'cif':

        with open(output_file_name, 'w') as f:
            simtk.openmm.app.pdbxfile.PDBxFile.writeFile(fixer.topology, fixer.positions, f, keepIds= True)

    elif output_file_type == 'pdb':

        with open(output_file_name, 'w') as f:
            simtk.openmm.app.PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds = True)
         
    return output_file_name
