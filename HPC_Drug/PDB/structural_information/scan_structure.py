######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################


#This file contains the functions necessary to scan the structure of a PDB file
#or a headerless mmCIF file

import copy


from HPC_Drug import important_lists
from HPC_Drug.auxiliary_functions import get_iterable
from HPC_Drug import orient


def get_metal_binding_residues_with_no_header(structure,
                                            cutoff = 3.0,
                                            protein_chain = 'A',
                                            protein_model = 0,
                                            COM_distance = 10.0,
                                            metals = important_lists.metals):

    """
    This function gets called by get_metalbinding_disulf_ligands

    This function iterates through the structure many times in order
    to return the metal binding residues through a substitution dictionary

    {residue_id : [residue_name, binding_atom, binding_metal]}

    It uses biopython structures

    structure :: a biopython structure of the protein

    cutoff :: double the maximum distance that a residue's center of mass and a metal ion
    can have to be considered binding default 3.0 angstrom

    protein_chain :: string default 'A', if == None no chain selection will be done

    protein_model :: integer default 0, if == None no model and no chain selection will be done

    metals :: a list (or tuple etc) that contains all the resnames (in capital letters) of metals necessary to look for,
    default HPC_Drug.important_lists.metals (Actually the easiest way to personalize metals is to append your custom values to this list)


    this function is slow and error prone
    and should only be used if there is no mmCIF with a good header
    
    It should not be necessary to change COM_distance because it simply is the distance between the center of mass of
    a residue and the metal that is used to know which atom distances to calculate
    """

        
    orient_object = orient.Orient()

    substitutions_dict = {}

    #selects model and chain if required
    if protein_model != None:

        try:

            model = structure[protein_model]

        except KeyError:

            model = structure

        # select only the right chain
        if protein_chain != None:

            try:

                chain = model[protein_chain.strip().upper()]

            except KeyError:

                chain = model

        else:
            chain = model

    else:
        chain = structure

    _chain = copy.deepcopy(chain)

    for residue in _chain:
        
        if residue.resname.strip().upper() in metals:
            for atom in residue:
                    
                #I get a second copy of all the residues in the chain
                tmp_struct = copy.deepcopy(chain)
                all_residues = tmp_struct.get_residues()

                #and iterate though them
                for other_residue in all_residues:
                    
                    #I avoid scanning the metal against it's self and against trash residues
                    if other_residue.resname.strip().upper() not in important_lists.metals:
                        COM_1, COM_2, distance = orient_object.center_mass_distance(structure_1 = residue, structure_2 = other_residue)
                        
                        if distance <= COM_distance:
                            
                            TMP_atom_dist = [1.E+20, 'DUMMY']
                            #check for the nearest atom of the binding residue
                            for other_atom in other_residue:

                                d = (atom.coord[0] - other_atom.coord[0])**2. + (atom.coord[1] - other_atom.coord[1])**2. + (atom.coord[2] - other_atom.coord[2])**2.
                                d = d ** (0.5)
                                
                                if d < TMP_atom_dist[0]:
                                    try:
                                        TMP_atom_dist = [d, other_atom.name.upper()]
                                    except:
                                        TMP_atom_dist = [d, other_atom.element.upper()]
                            
                            #checking if the nearest atom is near enough to be part of a binding residue
                            if TMP_atom_dist[0] <= cutoff:
                                #I add the other residue _id to the dictionary keys and give a value
                                substitutions_dict[other_residue.id[1]] = [other_residue.resname.strip().upper(), TMP_atom_dist[1], residue.resname.strip().upper()]

        
    #useless variables
    COM_1 = None
    COM_2 = None

    return substitutions_dict


def get_disulf_bonds_with_no_header(structure,
                                    cutoff = 3.0,
                                    protein_chain = 'A',
                                    protein_model = 0):
    """
    This function gets called by get_metalbinding_disulf_ligands
    
    This function iterates through the structure many times in order
    to return the disulf bonds through a substitution dictionary and a list of the binded couples

    {residue_id : [residue_name, binding_atom, binding_metal]}
    and
    [(cys_id, cys_id), (cys_id, ...), ...]

    return substitutions_dict, sulf_bonds

    it uses a biopython structure

    structure :: biopython structure of the protein

    cutoff :: double the maximum distance that two CYS S atoms
    can have to be considered binding default 3.0 angstrom

    protein_chain :: string default 'A', if == None no chain selection will be done

    protein_model :: integer default 0, if == None no model and no chain selection will be done

    this function is slow and error prone
    and should only be used if there is no mmCIF with a good header"""

    substitutions_dict = {}
    sulf_bonds = []

    #selects model and chain if required
    if protein_model != None:

        try:

            model = structure[protein_model]

        except KeyError:

            model = structure

        # select only the right chain
        if protein_chain != None:

            try:

                chain = model[protein_chain.strip().upper()]

            except KeyError:

                chain = model

        else:
            chain = model

    else:
        chain = structure

    _chain = copy.deepcopy(chain)

    for residue in _chain:
        
        if residue.resname.strip().upper() in important_lists.cyst_resnames:

            #I get a second copy of all the residues in the chain
            all_residues = copy.deepcopy(chain).get_residues()
            #and iterate though them
            for other_residue in all_residues:
                
                if other_residue.resname.strip().upper() in important_lists.cyst_resnames:
                    #I avoid scanning CYS against its self and agaisnt the already scanned ones
                    if other_residue.id[1] > residue.id[1]:

                        distance = (residue['SG'].coord[0] - other_residue['SG'].coord[0])**2. + (residue['SG'].coord[1] - other_residue['SG'].coord[1])**2. + (residue['SG'].coord[2] - other_residue['SG'].coord[2])**2.
                        distance = distance ** (0.5)
                        

                        if distance <= cutoff:

                            substitutions_dict[residue.id[1]] = ['CYS', 'SG', 'disulf']
                            substitutions_dict[other_residue.id[1]] = ['CYS', 'SG', 'disulf']

                            tmp_sulf_bond = [residue.id[1], other_residue.id[1]]

                            #In this way I avoid redundancy like [(1,2),(2,1)]
                            tmp_sulf_bond.sort()
                            tmp_sulf_bond = tuple(tmp_sulf_bond)

                            if tmp_sulf_bond not in sulf_bonds:

                                sulf_bonds.append(tmp_sulf_bond)

    return substitutions_dict, sulf_bonds

def get_organic_ligands_with_no_header(structure,
                                    protein_chain = 'A',
                                    protein_model = 0,
                                    trash = important_lists.trash,
                                    metals = important_lists.metals):
    """
    This function gets called by get_metalbinding_disulf_ligands

    This function iterates through the structure to get the organic ligand
    returning a list of lists containing
    [[resname, resnumber], [resname, resnumber], ...]
    If there are none returns None

    it uses a biopython structure

    structure :: biopython structure of the protein

    protein_chain :: string default 'A', if == None no chain selection will be done

    protein_model :: integer default 0, if == None no model and no chain selection will be done

    trash :: a list (or tuple etc) that contains all the resnames (in capital letters) of trash ligands to avoid listing,
    default HPC_Drug.important_lists.trash (Actually the easiest way to personalize trash is to append your custom values to this list)

    metals :: a list (or tuple etc) that contains all the resnames (in capital letters) of metals necessary to look for,
    default HPC_Drug.important_lists.metals (Actually the easiest way to personalize metals is to append your custom values to this list)

    this function is slow and error prone
    and should only be used if there is no mmCIF with a good header"""

    ligand_list = []

    #selects model and chain if required
    if protein_model != None:

        try:

            model = structure[protein_model]

        except KeyError:

            model = structure

        # select only the right chain
        if protein_chain != None:

            try:

                chain = model[protein_chain.strip().upper()]

            except KeyError:

                chain = model

        else:
            chain = model

    else:
        chain = structure


    for residue in chain:

        #it's an hetero atom
        if residue.id[0].strip() != '':

            if residue.resname.strip().upper() not in metals and residue.resname.strip().upper() not in trash:

                ligand_list.append([residue.resname.strip().upper(), residue.id[1]])

    #if there are no ligand returns None
    if len(ligand_list) == 0:

        print("The list of ligands is empty, going on returning a None item")
        ligand_list = None

    return ligand_list


def get_metalbinding_disulf_ligands(Protein, trash = important_lists.trash, metals = important_lists.metals):
    """
    This is a template that uses the other funcions on this file to return a dictionary
    with key = resnum and value = (resname, binding atom, metal) or ('CYS', 'SG', 'disulf')
    depending if the residue number {resnum} binds a metallic ion or is part of a disulf bond
    and updates Protein.substitutions_dict with it

    and a list of tuples that contain the couples of CYS that are part of a disulf bond
    and updates Protein.sulf_bonds with it

    and a list of tuples with the residue name and residue number of the organic ligands
    (if there are none None will be returned)

    Protein :: a HPC_Drug.structures.protein.Protein instance

    trash :: a list (or tuple etc) that contains all the resnames (in capital letters) of trash ligands to avoid listing,
    default HPC_Drug.important_lists.trash (Actually the easiest way to personalize trash is to append your custom values to this list)

    metals :: a list (or tuple etc) that contains all the resnames (in capital letters) of metals necessary to look for,
    default HPC_Drug.important_lists.metals (Actually the easiest way to personalize metals is to append your custom values to this list)

    return Protein, [[lig_resname, lig_resnum], ...]
    """

    Protein.update_structure(struct_type = "biopython")

    #metal binding residues
    metal_binding_dict = get_metal_binding_residues_with_no_header(structure = Protein.structure,
                                                                protein_chain = Protein.chain,
                                                                protein_model = Protein.model,
                                                                metals = metals)

    #disulf bonds
    disulf_binding_dict, Protein.sulf_bonds = get_disulf_bonds_with_no_header(structure = Protein.structure,
                                                                            protein_chain = Protein.chain,
                                                                            protein_model = Protein.model)


    #I don't want the disulf_binding_dict to overwrite the metal_binding_dict elements during the merging
    #because of course if there are multiple CYS binding one metal they will be near enough
    #to be included in the disulf_binding_dict (false positives)
    for key in disulf_binding_dict.keys():

        if key not in metal_binding_dict.keys():

            metal_binding_dict[key] = disulf_binding_dict[key]


    Protein.substitutions_dict = metal_binding_dict

    organic_ligand_list = get_organic_ligands_with_no_header(structure = Protein.structure,
                                                            protein_chain = Protein.chain,
                                                            protein_model = Protein.model,
                                                            trash = trash,
                                                            metals = metals)

    return Protein, organic_ligand_list