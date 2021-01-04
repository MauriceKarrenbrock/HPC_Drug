######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

import Bio.PDB

def remove_disordered_atoms(Protein):
    """
    Removes disordered atoms, solves a problem about "copied atoms don't inherit disordered_get_list in Biopython"

    Protein :: HPC_Drug.structures.protein.Protein instance

    return Protein
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

    if Protein.file_type == 'pdb':
        p = Bio.PDB.PDBParser()

    elif Protein.file_type == 'cif':
        p = Bio.PDB.MMCIFParser()

    else:
        raise TypeError(f"Protein.file_type must be of 'cif' or 'pdb' type not {Protein.file_type}")
   
    Protein.structure = p.get_structure(Protein.protein_id, Protein.pdb_file)

    if Protein.file_type == 'pdb':
        s = Bio.PDB.PDBIO()

    elif Protein.file_type == 'cif':
        s = Bio.PDB.MMCIFIO()

    s.set_structure(Protein.structure)
    s.save(Protein.pdb_file)

    return Protein