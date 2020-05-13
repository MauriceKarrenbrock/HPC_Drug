######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
Remember to use the orac residue renaming tool before you use this functions
"""

from HPC_Drug.MD.orac import get_resnum_cutoff
from HPC_Drug import important_lists

def custom_orac_seqres_from_PDB(Protein):

    """
    Updates Protein.seqres with the seqres format needed for the
    Orac input

    Protein :: HPC_Drug.structures.structure.Protein instance
    
    returns Protein
    """

    Protein.update_structure("biopython")

    Protein.seqres = []
    tmp_metals = []

    for model in Protein.structure:
        for chain in model:

            #get residue list of every chain
            tmp_seqres = []

            for residue in chain:

                #skipping HETATMs
                if residue.id[0].strip() == '' and residue.resname.strip() not in important_lists.metals:
                    
                    #In the pdb this residues are called CYM because Biopython structures can
                    #only have 3 letters resnames, but Orac's tpg file calls it CYSM
                    #So I am modifying it in the seqres list
                    #It's probably not the best way to do it
                    if residue.resname.strip() == 'CYM':
                        tmp_seqres.append('CYSM')
                    else:
                        tmp_seqres.append(residue.resname.strip())
                
                #putting metal ions and hetatms in a separated list
                else:
                    tmp_metals.append(residue.resname.strip())

            
            #rename first and last residue of every chain
            tmp_seqres[0] = tmp_seqres[0] + '-H'
            tmp_seqres[-1] = tmp_seqres[-1] + '-O'

            for tmp in tmp_seqres:
                Protein.seqres.append(tmp)

        #The residue id of the first residue (it is not always one)
        cutoff = get_resnum_cutoff.get_resnum_cutoff(structure = Protein.structure)

        for i, resname in enumerate(Protein.seqres):

            if resname == 'CYS':
                if not (str(i + cutoff) in Protein.substitutions_dict.keys()):

                    Protein.seqres[i] = 'CYSH'
    
    for tmp in tmp_metals:
        Protein.seqres.append(tmp)
    
    return Protein
