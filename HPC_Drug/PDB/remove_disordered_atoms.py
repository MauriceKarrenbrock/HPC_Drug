######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

import Bio.PDB

def remove_disordered_atoms(input_pdb, output_pdb):
    """
    Removes disordered atoms, solves a problem about "copied atoms don't inherit disordered_get_list in Biopython"

    Parameters
    ------------
    input_pdb : str
        a pdb or mmcif file
    output_pdb : str
        a pdb or mmcif file

    return
    ----------
    output_pdb : str
    """

    # This part eliminates any disordered atom part,
    #Only keeps one possible position
    #because copied atoms don't inherit disordered_get_list in Biopython
    def get_unpacked_list(self):
        """
        Returns all atoms from the residue,
        in case of disordered, keep only first alt loc and remove the alt-loc tag
        """

        atom_list = self.get_list()
        undisordered_atom_list = []
        for atom in atom_list:
            if atom.is_disordered():
                atom.altloc=" "
                undisordered_atom_list.append(atom)
            else:
                undisordered_atom_list.append(atom)
        return undisordered_atom_list


    Bio.PDB.Residue.Residue.get_unpacked_list = get_unpacked_list

    if str(input_pdb).split('.')[-1] == 'pdb':
        p = Bio.PDB.PDBParser()

    elif str(input_pdb).split('.')[-1] == 'cif':
        p = Bio.PDB.MMCIFParser()

    else:
        raise TypeError(f"input_pdb must be of 'cif' or 'pdb' type not {str(input_pdb).split('.')[-1]}")
   
    structure = p.get_structure('aaaa', str(input_pdb))

    if str(output_pdb).split('.')[-1] == 'pdb':
        s = Bio.PDB.PDBIO()

    elif str(output_pdb).split('.')[-1] == 'cif':
        s = Bio.PDB.MMCIFIO()
    
    else:
        raise TypeError(f"output_pdb must be of 'cif' or 'pdb' type not {str(input_pdb).split('.')[-1]}")

    s.set_structure(structure)
    s.save(str(output_pdb))

    return output_pdb