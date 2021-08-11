######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

import openbabel

def convert_and_protonate_pdb_to_sdf(pdb_file, sdf_file, ph=7.0, ligand_resname='LIG'):
    """Uses openbabel to convert a pdb to ad sdf and adds hydrogens at the given pH

    Parameters
    ------------
    pdb_file : str or path
        the input pdb file
    sdf_file : str or path
        the output sdf file
    ph : float, default=7.0
        if set to None hydrogens won't be touched
    ligand_resname : str, default='LIG'

    Notes
    ----------
    when starting from a pdb it is better to use openbabel=2.4
    and not openbabel=3.x
    see:
    https://sourceforge.net/p/openbabel/mailman/openbabel-discuss/thread/CAOC-GK0-98HUDpGSRkyZ8G4pKoXyawNK8ZrvVtJLEDwYBV19SQ%40mail.gmail.com/
    """

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "sdf")

    mol = openbabel.OBMol()

    obConversion.ReadFile(mol, str(pdb_file))

    if ph is not None:
        mol.DeleteHydrogens()
        mol.AddHydrogens(False, True, ph)

    obConversion.WriteFile(mol, str(sdf_file))

    # open babel gives some random names
    # so I change the first line
    with open(sdf_file, 'r') as f:
        lines = f.readlines()
        
    lines[0] = ligand_resname.strip() + '\n'
    with open(sdf_file, 'w') as f:
        for line in lines:
            f.write(line)
