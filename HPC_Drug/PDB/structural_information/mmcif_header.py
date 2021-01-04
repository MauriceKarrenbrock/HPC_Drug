######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################


#This file contains the functions necessary to parse the header of a mmCIF file

from HPC_Drug import important_lists
from HPC_Drug.auxiliary_functions import get_iterable

def get_metal_binding_residues(mmcif2dict, metals = important_lists.metals):
    """
    This function is called from get_metalbinding_disulf_ligands

    Searces the given mmcif file for the metal binding residues parsing the header
    returns a dictionary that has as key the residue number and as vaue a tuple with
    (resname, binding atom, metal) for metal binding residues

    mmcif2dict :: a dictionary of the type you obtain with HPC_Drug.PDB.biopython.mmcif2dict function

    metals :: a list (or tuple etc) that contains all the resnames (in capital letters) of metals necessary to look for,
    default HPC_Drug.important_lists.metals (Actually the easiest way to personalize metals is to append your custom values to this list)

    return {resnum : (resname, binding atom, metal), ...}
    """

    substitutions = {}

    #working variable
    _mmcif2dict = {}

    #trasforms all the not iterable values of the dictionary in
    #iterable tuples
    for key in mmcif2dict.keys():
        _mmcif2dict[key] = get_iterable.get_iterable(mmcif2dict[key])

    #I check for residues binding metals
    #In the end I get a dictionary that has as key the residue number
    #and as vaue a tuple with (resname, binding atom, metal) for metal binding residues
    if '_struct_conn.conn_type_id' in _mmcif2dict.keys():
        for i, bound_type in enumerate(_mmcif2dict['_struct_conn.conn_type_id']):

            if bound_type == 'metalc':

                if _mmcif2dict['_struct_conn.ptnr1_label_comp_id'][i]\
                    in metals\
                    or _mmcif2dict['_struct_conn.ptnr2_label_comp_id'][i]\
                    in metals:

                    if _mmcif2dict['_struct_conn.ptnr1_label_comp_id'][i] in metals:
                        metal_index = '1'
                        res_index = '2'
                    elif _mmcif2dict['_struct_conn.ptnr2_label_comp_id'][i] in metals:
                        metal_index = '2'
                        res_index = '1'

                    # substitutions[resnum] = (res_label_comp_id, res_label_atom_id, metal_label_comp_id)
                    substitutions[int(_mmcif2dict[f'_struct_conn.ptnr{res_index}_auth_seq_id'][i].strip())] = (
                                    _mmcif2dict[f'_struct_conn.ptnr{res_index}_label_comp_id'][i],
                                    _mmcif2dict[f'_struct_conn.ptnr{res_index}_label_atom_id'][i],
                                    _mmcif2dict[f'_struct_conn.ptnr{metal_index}_label_comp_id'][i]
                                    )
                
                else:
                    string = f"Not implemented metal bound, going on pretending nothing happened\n\
                    More info:\n\
                    _struct_conn.ptnr1_label_comp_id = \
                    {_mmcif2dict['_structure_conn.ptnr1_label_comp_id'][i]}\n\
                    _struct_conn.ptnr2_label_comp_id = \
                    {_mmcif2dict['_structure_conn.ptnr2_label_comp_id'][i]}" 

                    print(string)

            else:
                pass

    return substitutions


def get_disulf_bonds(mmcif2dict):
    """
    This function is called from get_metalbinding_disulf_ligands

    Searces the given mmcif file for disulf bonds parsing the header
    returns a dictionary that has as key the residue number and as vaue a tuple with
    ('CYS', 'SG', 'disulf') for any disulf cysteine.
    And a a list composed of tuples containing the resnumbers of the 2 CYS that bound through disulfide bond

    mmcif2dict :: a dictionary of the type you obtain with HPC_Drug.PDB.biopython.mmcif2dict function

    return {resnum : ('CYS', 'SG', 'disulf'), ...}   [(resnum, resnum), (...), ...]
    """

    sulf_bonds = []
    substitutions = {}

    #working variable
    _mmcif2dict = {}

    for key in mmcif2dict.keys():
        _mmcif2dict[key] = get_iterable.get_iterable(mmcif2dict[key])

    #Create a tuple with ('CYS', 'SG', 'disulf') for the disulfide bonds.
    #Separately I create a list composed of tuples containing the resnumbers
    #of the 2 CYS that bound through disulfide bond
    if '_struct_conn.conn_type_id' in _mmcif2dict.keys():
        for i, bound_type in enumerate(_mmcif2dict['_struct_conn.conn_type_id']):

            if bound_type == 'disulf':

                #if it is present ptnr1_auth_seq_id is the resnumber that remains when transforming it in a pdb
                try:

                    sulf_bond_tmp = [int(_mmcif2dict['_struct_conn.ptnr1_auth_seq_id'][i].strip()),
                                    int(_mmcif2dict['_struct_conn.ptnr2_auth_seq_id'][i].strip())]

                except:

                    sulf_bond_tmp = [int(_mmcif2dict['_struct_conn.ptnr1_label_seq_id'][i].strip()),
                                    int(_mmcif2dict['_struct_conn.ptnr2_label_seq_id'][i].strip())]

                sulf_bond_tmp.sort()

                sulf_bond_tmp = tuple(sulf_bond_tmp)

                #removes redundancy of type : [(1,2),(2,1)]
                if sulf_bond_tmp not in sulf_bonds:

                    sulf_bonds.append(sulf_bond_tmp)
                
                substitutions[sulf_bond_tmp[0]] = ('CYS', 'SG', 'disulf')
                substitutions[sulf_bond_tmp[1]] = ('CYS', 'SG', 'disulf')

            else:
                pass

    return substitutions, sulf_bonds


def get_organic_ligands(mmcif2dict, protein_chain = None, trash = important_lists.trash, metals = important_lists.metals):
    """
    This function is called from get_metalbinding_disulf_ligands

    Searces the given mmcif file for organic ligands parsing the header
    returns a list of resnames. If protein chain is None (default) will list all ligands from
    any chain, if protein chain is set does only consider the ones of the given chain (es A)

    mmcif2dict :: a dictionary of the type you obtain with HPC_Drug.PDB.biopython.mmcif2dict function

    protein_chain :: string, default None, the chain id of the chain you want to analize in capital letters (es A)

    trash :: a list (or tuple etc) that contains all the resnames (in capital letters) of trash ligands to avoid listing,
    default HPC_Drug.important_lists.trash (Actually the easiest way to personalize trash is to append your custom values to this list)

    metals :: a list (or tuple etc) that contains all the resnames (in capital letters) of metals necessary to look for,
    default HPC_Drug.important_lists.metals (Actually the easiest way to personalize metals is to append your custom values to this list)

    return [resname, resname, ...]
    """

    if protein_chain != None:
        protein_chain = protein_chain.strip().upper()

    ligand_resnames = []

    if '_struct_site.details' in mmcif2dict.keys():

        #record the ligands resname
        for i in mmcif2dict['_struct_site.details']:
            n = i.split()

            #if required
            #chooses the right chain only
            if protein_chain != None:
                dont_check_chain = False
            else:
                dont_check_chain = True

            if n[-2].strip() == protein_chain or dont_check_chain:
                n = n[-3]
                n = n.strip()

                if n not in metals:
                    if n not in trash:
                        ligand_resnames.append(n)

        #Removes redundancy
        ligand_resnames = set(ligand_resnames)
        ligand_resnames = list(ligand_resnames)

    return ligand_resnames


def get_ligand_resnum(structure, ligand_resnames = None, protein_chain = 'A', protein_model = 0):
    """
    This function is called from get_metalbinding_disulf_ligands

    Given a Biopython structure and a list of Ligand_resnames will return 
    a list containing the ligand resnames and resnumbers
    in order to distinguish ligands with the same resname: [[resname, resnumber], [.., ...], ...]
    
    ligand_resnames :: list, it is a list containing the organic ligand resnames (capital letters) to look for
    if it is == None or empty will return None

    protein_chain :: string, default A, the chain id of the chain you want to analize in capital letters (es A), if == None no chain selection will be done
    
    protein_model :: integer, default 0, the model to check, if == None no chain and no model selection will be done

    return [[resname, resnumber], [.., ...], ...]
    """
    
    # If the ligand resname is a single string I transform it in an iterable object
    ligand_resnames = get_iterable.get_iterable(ligand_resnames)

    if len(ligand_resnames) == 0 or ligand_resnames[0] == None:
        print("The list of ligands is empty, going on returning a None item")
        return None

    if protein_model != None:

        try:

            _structure = structure[protein_model]

        except KeyError:
                _structure = structure

        if protein_chain != None:
            try:
                #Taking only the right chain and model
                _structure = _structure[protein_chain]
            except KeyError:
                _structure = structure
    
    else:
        _structure = structure

    residues = _structure.get_residues()

    ligand_list = []

    for residue in residues:
        if (residue.resname.strip().upper() in ligand_resnames) or (residue.resname.strip().lower() in ligand_resnames):
            ligand_list.append([residue.resname.strip(), residue.id[1]])

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

    Protein.update_structure(struct_type = "mmcif2dict")

    #metal binding residues
    metal_binding_dict = get_metal_binding_residues(mmcif2dict = Protein.structure, metals = metals)

    #disulf bonds
    disulf_binding_dict, Protein.sulf_bonds = get_disulf_bonds(mmcif2dict = Protein.structure)

    Protein.substitutions_dict = {**metal_binding_dict, **disulf_binding_dict}

    ligand_resnames = get_organic_ligands(mmcif2dict = Protein.structure, protein_chain = Protein.chain, metals = metals, trash = trash)

    Protein.update_structure(struct_type = "biopython")

    organic_ligand_list = get_ligand_resnum(
                                    structure = Protein.structure,
                                    ligand_resnames = ligand_resnames,
                                    protein_chain = Protein.chain,
                                    protein_model = Protein.model)

    return Protein, organic_ligand_list



