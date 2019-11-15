#containing functions and classes for the creation of a gromacs input

import structures
import Bio.PDB
import file_manipulation
import important_lists

def residue_substitution(Protein, substitution = 'standard'):
    """Takes a protein instance, and returns one

    renames the resnames of the metal binding residues in order
    to generate the right force field"""

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

    p = Bio.PDB.PDBParser()
    struct = p.get_structure(Protein.protein_id, Protein.filename)

    #I iterate through the structure
    for model in struct:
        for chain in model:
            for residue in chain:

                #search for the right residues
                res_id = str(residue._id[1])

                #check if they are bounding a metal and are not a disulfide bond
                if res_id in Protein.substitutions_dict.keys():
                    if Protein.substitutions_dict[res_id][2] in important_lists.metals:

                        #make substitutions with the selected function
                        residue.resname = substitute(Protein.substitutions_dict[res_id])

    Protein.structure = struct
    Protein.filename = Protein.write_PDB(Protein.filename, 'biopython')

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

    elif resname == 'CIM' and input_list[2] == 'ZN':
        resname = 'CIZ'

    return resname