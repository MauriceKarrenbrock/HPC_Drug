######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
Residue renaming for gromacs
"""

from HPC_Drug import important_lists


def residue_substitution(Protein, substitution = 'standard', ph = 7.0):
    """Takes a protein instance, and returns one

    renames the resnames of the metal binding residues in order
    to generate the right force field
    
    decides how histidines should be protonated (ph > 5 HSD, ph < 5 HIS except for the metal binding ones)"""

    #You can choose any implemented substitution
    #or create your own
    if substitution == 'standard':
        def substitute(input_list) :
            return standard_gromacs_substitutions(input_list)
    
    elif substitution == 'custom_zinc':
        def substitute(input_list) :
            return custom_zinc_gromacs_substitutions(input_list)
    
    else:
        raise NotImplementedError(substitution)

    Protein.update_structure(struct_type = "biopython")


    #I iterate through the structure
    for model in Protein.structure:
        for chain in model:
            for residue in chain:

                #search for the right residues
                res_id = residue.id[1]

                #check if they are bounding a metal and are not a disulfide bond
                if res_id in Protein.substitutions_dict.keys():
                    if Protein.substitutions_dict[res_id][2] in important_lists.metals:

                        #make substitutions with the selected function
                        residue.resname = substitute(Protein.substitutions_dict[res_id])

                    #if a cysteine is in a disulfide bond it is called CYX
                    elif residue.resname in important_lists.cyst_resnames:
                        residue.resname = _disulf_cysteine_gromacs_substitutions(
                                                    input_list = Protein.substitutions_dict[res_id],
                                                    resname = residue.resname
                                                    )

                #give the right name to histidines in order to get the right protonation
                elif residue.resname in important_lists.hist_resnames:
                    
                    if ph < 6.0:
                        residue.resname = 'HIS'
                    
                    else:
                        residue.resname = 'HID'


    Protein.write(file_name = Protein.pdb_file, struct_type = 'biopython')

    return Protein


def standard_gromacs_substitutions(input_list):
    """It is called by residue_substitution if run with standard
    substitution option
    takes a tuple long 3 with (resname, binding atom, metal name)
    
    returns the right resname"""
    
    if input_list[0] in important_lists.hist_resnames:

        if input_list[1] == 'NE2':
            resname = 'HID'
        
        elif input_list[1] == 'ND1':
            resname = 'HIE'
        
        else:
            raise NotImplementedError(f"{input_list[1]} is not an implemented atom name for histidine")

    elif input_list[0] in important_lists.cyst_resnames:
        resname = 'CYM'

    else:
        print(input_list[0], 'substitution not implemented, going on like nothing happened')
        resname = input_list[0]

    return resname

def custom_zinc_gromacs_substitutions(input_list):
    """It is a custom version of standard_gromacs_substitutions(list) function
    with custom names for zinc binding residues"""
    
    resname = standard_gromacs_substitutions(input_list)

    if resname == 'HID' and input_list[2] == 'ZN':
        resname = 'HDZ'
    
    elif resname == 'HIE' and input_list[2] == 'ZN':
        resname = 'HEZ'

    elif resname == 'CYM' and input_list[2] == 'ZN':
        resname = 'CYZ'

    return resname

def _disulf_cysteine_gromacs_substitutions(input_list, resname = None):
    """Gives the CYX resname to the disulfide bond cysteines"""
    
    if input_list[2].lower().strip() == "disulf":

        return 'CYX'

    elif resname != None:

        return resname

    else:

        return input_list[0]
