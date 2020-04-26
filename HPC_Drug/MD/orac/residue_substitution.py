######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
The functions to rename the residues accordingly to Orac
"""

from HPC_Drug import important_lists

def residue_substitution(Protein, substitution = 'standard', ph = 7.0):
    """
    Takes a HPC_Drug.structures.protein.Protein instance, and returns one

    renames the resnames of the metal binding residues
    and of the cysteines in order
    to generate the right force field
    
    decides how histidines should be protonated (ph > 5 HSD, ph < 5 HIS except for the metal binding ones)
    
    Protein :: HPC_Drug.structures.protein.Protein

    substitution :: string, default "standard", options "standard" "custom_zinc" 
    it is the kind of renaming you are interested in

    ph :: float, dafault 7.0
    """

    #You can choose any implemented substitution
    #or create your own
    if substitution == 'standard':
        def substitute(input_list) :
            return standard_orac_substitutions(input_list)
    
    elif substitution == 'custom_zinc':
        def substitute(input_list) :
            return custom_zinc_orac_substitutions(input_list)
    
    else:
        raise NotImplementedError(substitution)
    
    Protein.update_structure(struct_type = "biopython")

    #I iterate through the structure
    for model in Protein.structure:
        for chain in model:
            for residue in chain:

                #search for the right residues
                res_id = str(residue.id[1])

                #check if they are bounding a metal and are not a disulfide bond
                if res_id in Protein.substitutions_dict.keys():
                    if Protein.substitutions_dict[res_id][2] in important_lists.metals:

                        #make substitutions with the selected function
                        residue.resname = substitute(Protein.substitutions_dict[res_id])

                #give the right name to histidines in order to get the right protonation
                elif residue.resname in important_lists.hist_resnames:
                    
                    if ph < 6.0:
                        residue.resname = 'HIS'
                    
                    else:
                        residue.resname = 'HSD'

    Protein.write(file_name = Protein.pdb_file, struct_type = 'biopython')

    return Protein

def standard_orac_substitutions(input_list):
    """
    It is called by residue_substitution if run with standard
    substitution option
    takes a tuple long 3 with (resname, binding atom, metal name)
    
    returns the right resname
    """
    
    if input_list[0] in important_lists.hist_resnames:

        if input_list[1] == 'NE2':
            resname = 'HSD'
        
        elif input_list[1] == 'ND1':
            resname = 'HSE'
        
        else:
            raise NotImplementedError(f"{input_list[1]} is not an implemented atom name for histidine")

    elif input_list[0] in important_lists.cyst_resnames:
        #This isn't the right name to write in the input seqres for Orac
        #But as Biopython strucures can only contain 3 letter resnames
        #I am giving this as a temporary name
        #It will be changed to CYSM in the get seqres function
        #This part shall be redone better
        resname = 'CYM'

    else:
        print(input_list[0], 'substitution not implemented, going on like nothing happened')
        resname = input_list[0]

    return resname

def custom_zinc_orac_substitutions(input_list):
    """
    Called by residue_substitution

    It is a custom version of standard_gromacs_substitutions(list) function
    with custom names for zinc binding residues
    """
    
    resname = standard_orac_substitutions(input_list)

    if resname == 'HSD' and input_list[2] == 'ZN':
        resname = 'HDZ'
    
    elif resname == 'HSE' and input_list[2] == 'ZN':
        resname = 'HEZ'

    elif resname == 'CYM' and input_list[2] == 'ZN':
        resname = 'CYZ'

    return resname
