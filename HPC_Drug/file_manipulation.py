######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################


"""
DEPRECATED is here only because old stuff still uses it

This file contains the classes for file manipulation like PDB & mmCIF
"""

from HPC_Drug import important_lists
from HPC_Drug import orient
from HPC_Drug.auxiliary_functions import get_iterable


import Bio.PDB
import Bio.PDB.MMCIF2Dict
import prody


def mmcif2pdb(Protein = None):
    """Takes a Protein instance and returns one
    if Protein.file_type == 'cif' converts the file in a pdb
    and updates Protein.file_type and Protein.pdb_file
    otherwise does nothing"""

    if Protein == None:
        raise TypeError("I need a PRotein as input, not a None type")

    if Protein.file_type == 'pdb':
        return Protein

    elif Protein.file_type == 'cif':

        new_name = Protein.pdb_file.rsplit('.', 1)[0] + '.pdb'

        p = Bio.PDB.MMCIFParser()

        struct = p.get_structure(Protein.protein_id, Protein.pdb_file)

        s = Bio.PDB.PDBIO()
        s.set_structure(struct)
        s.save(new_name)

        Protein.pdb_file = new_name
        Protein.file_type = 'pdb'

        return Protein

    else:
        raise ValueError(f"Protein.file_type must be 'pdb' or 'cif' not {Protein.file_type}")


class FileCruncer(object):
    """Generic file cruncer class"""
    def __init__(self):
        raise NotImplementedError('This kind file was not implemented yet')

def ProteinCruncer(file_type):
    """Returns the right file cruncer class
    for the given file_type"""

    if file_type == 'pdb':
        return PDBCruncer()
    elif file_type == 'cif':
        return MMCIFCruncer()
    else:
        raise NotImplementedError(f"This kind file was not implemented yet: {file_type}, only 'pdb' and 'cif'")


class PDBCruncer(FileCruncer):
    """Contains the metods for working on PDB
    with the ProDy module"""

    def __init__(self):
        pass

    def parse(self, protein_id = None, filename = None):
        """Parses a PDB with the ProDy parser"""
        
        
        if filename == None:

            #check that what I got is a string
            if type(protein_id) != str:
                raise TypeError(f"Need a string type, not a {type(protein_id)}")

            parser = prody.parsePDB(protein_id)
        else:
            #check that what I got is a string
            if type(filename) != str:
                raise TypeError(f"Need a string type, not a {type(filename)}")

            parser = prody.parsePDB(filename)

        return parser
        
    def get_protein(self,
                    protein_id = None,
                    filename = None,
                    structure = None):
        
        """returns a protein ProDy structure"""


        if structure == None:
            parser = self.parse(protein_id, filename)
        else:
            parser = structure

        try:
            protein_structure = parser.select('protein or ion')
        except:
            protein_structure = parser.select('protein')

        return protein_structure

    def get_ligand(self,
                protein_id = None,
                filename = None,
                ligand_name = None,
                structure = None,
                ligand_resnumber = None):
        """Creates a PDB the organic ligand (if any)
        with ProDy
        If structure is given it must be a prody structure"""

        if structure == None:
            parser = self.parse(protein_id, filename)
        
        elif not isinstance(structure, prody.AtomGroup):
            
            raise TypeError(f"Need a Prody structure (prody.AtomGroup)\n not a {type(structure)}")

        else:
            parser = structure

        if ligand_resnumber != None and type(ligand_resnumber) != str:
            ligand_resnumber = str(ligand_resnumber)

        #if no resnumber is given the selection will be made using the resname
        #It makes impossible to distinguish between ligands with the same resname
        if ligand_resnumber == None and ligand_name != None and len(ligand_name) != 0:
            ligand_structure = parser.select('resname ' + ligand_name)
            print(f"Selected ligand {ligand_name}")

        elif ligand_resnumber != None:
            ligand_structure = parser.select('resnum ' + ligand_resnumber)
            print(f"Selected ligand resnumber: {ligand_resnumber}")
        
        else:
            print("No ligand found will return None")
            return None
        
        return ligand_structure

    
    def write_PDB(self, structure, pdb_name):
        """
        Write a prody protein to a pdb file
        :param protein: protein object from prody
        :param pdb_name: base name for the pdb file
        :return: None
        """
        
        prody.writePDB(f"{pdb_name}", structure)



class MMCIFCruncer(FileCruncer):
    """Contains the metods for working on mmCIF
    with th ProDy module"""

    def __init__(self):
        pass

    def parse(self, protein_id = None, filename = None, Protein_model = None):
        """Parses a PDB with the ProDy parser"""
        
        if filename == None:
            parser = prody.parseCIF(protein_id, model = Protein_model)
        else:
            parser = prody.parseCIF(filename, model = Protein_model)

        return parser
    
    def get_protein(self, protein_id = None, filename = None, Protein_model = None, structure = None):

        if structure == None:
            parser = self.parse(protein_id, filename, Protein_model)
        else:
            parser = structure

        p = PDBCruncer()

        protein = p.get_protein(protein_id, filename, parser)

        return protein

    def get_ligand(self,
                protein_id = None,
                filename = None,
                ligand_name = None,
                structure = None):

        """Creates a PDB with all the not water HETATM
        with ProDy"""

        if structure == None:
            parser = self.parse(protein_id, filename)
        else:
            parser = structure

        p = PDBCruncer()

        ligand = p.get_ligand(protein_id, filename, parser)

        return ligand
    
    def write_MMCIF(self, structure, pdb_name):
        """
        Not implementes, writes a PDB instead
        using PDBCruncer() class
        """
        print("Not implemented, printing PDB instead")

        p = PDBCruncer()
        p.write_PDB(f"{pdb_name}", structure)

class PDBRepair(FileCruncer):
    """This class contains the various options to add
    missing atoms and missing residues,
    and change all non standard residues to standard ones"""

    def __init__(self):
        pass
    
    def add_missing_atoms(self, pdb_id, input_filename,
                        file_type = 'cif',
                        repairing_method = 'pdbfixer',
                        output_filename = None,
                        ph = 7.0,
                        add_H = False):
        """This function returns a repaired PDBx/mmCIF with missing atoms but not hydrogens
        chooses the right repairing method with repairing_method (pdbfixer is default)
        The returned PDB is called output_filename
        if output_filename is None (default) the file is saved as
        {protein_id}_repaired.cif in the working directory
        returns the output_filename

        does only work on the protein, not on the ligand
        
        The input file MUST be a PDBx/mmCIF or PDB"""

        if repairing_method == 'pdbfixer':
            import pdbfixer
            import simtk.openmm.app
            import os

            if file_type == 'cif':


                #PART OF A PATCH TO ADRESS issue #195 on pdbfixer github
                TMP_dict = Bio.PDB.MMCIF2Dict.MMCIF2Dict(input_filename)

                patch_dict = {}

                for i in range(len(TMP_dict['_atom_site.label_asym_id'])):

                    patch_dict[TMP_dict['_atom_site.label_asym_id'][i]] = TMP_dict['_atom_site.auth_asym_id'][i]
                #--------

            fixer = self.fix_pdbfixer(input_filename = input_filename,
                                    ph = ph,
                                    add_H = add_H,
                                    file_type = file_type)

            if output_filename == None:
                output_filename = f'{pdb_id}_repaired.{file_type}'
            
            if file_type == 'cif':

                with open(output_filename, 'w') as f:
                    simtk.openmm.app.pdbxfile.PDBxFile.writeFile(fixer.topology, fixer.positions, f, keepIds= True)


                #PART OF A PATCH TO ADRESS issue #195 on pdbfixer github
                TMP_dict = Bio.PDB.MMCIF2Dict.MMCIF2Dict(output_filename)

                for i in range(len(TMP_dict['_atom_site.label_asym_id'])):

                    TMP_dict['_atom_site.auth_asym_id'][i] = patch_dict[TMP_dict['_atom_site.label_asym_id'][i]]
                
                p = Bio.PDB.MMCIFIO()
                p.set_dict(TMP_dict)
                p.save(output_filename)
                #---------------


            elif file_type == 'pdb':

                with open(output_filename, 'w') as f:
                    simtk.openmm.app.PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds = True)

            else:
                raise ValueError(f"file_type must be cif or pdb not {file_type}")
            
            output_filename = os.getcwd() + '/' + output_filename
            if not os.path.exists(output_filename):
                raise Exception(f'was not able to write {output_filename}')
            else:
                
                
                return output_filename
        
        else:
            raise NotImplementedError(f'{repairing_method}')


    def fix_pdbfixer(self, input_filename, ph = 7.0, add_H = False, file_type = 'cif'):
        """This method is called by add_missing_atoms
        Reads an PDBx/mmCIF or PDB"""
        
        import pdbfixer
        import simtk.openmm.app

        if file_type == 'cif':
            with open(input_filename, 'r') as f:
                fixer = pdbfixer.PDBFixer(pdbxfile = f)

        elif file_type == 'pdb':

            with open(input_filename, 'r') as f:
                fixer = pdbfixer.PDBFixer(pdbfile = f)
        
        else:
            raise ValueError(f"file_type must be cif or pdb not {file_type}")

        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        #fixer.removeHeterogens(False)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        if add_H == True:
            fixer.addMissingHydrogens(ph)
        #fixer.addSolvent(fixer.topology.getUnitCellDimensions())

        return fixer


class SubstitutionParser(FileCruncer):

    def __init__(self):
        pass
    
    def parse_substitutions_PDB(self, file_name, protein_chain):

        """Checks the necessary substitutions to the residue names
        needed for the correct placing of the hydrogen atoms
        for a MD (by the MD program)
        (very important when there are coordinated metals)
        returns a vocabulary for substitutions
        key = residue_id
        value = new_name

        and if it's present the resname used for the organic ligand (list)
        """
        
        substitutions = {}
        sulf_bonds = []

        metals = important_lists.metals
        
        #list of trash ligands
        trash = important_lists.trash

        cif_dict = Bio.PDB.MMCIF2Dict.MMCIF2Dict(file_name)

        #trasforms all the not iterable values of the dictionary in
        #iterable tuples
        for key in cif_dict.keys():
            cif_dict[key] = get_iterable.get_iterable(cif_dict[key])

        ligand = self.parse_ligands_from_header(cif_dict = cif_dict, metals = metals, trash = trash, protein_chain = protein_chain)         
    

        #I check for residues binding metals and disulfide bonds
        #In the end I get a dictionary that has as key the residue number
        #and as vaue a tuple with (resname, binding atom, metal) for metal binding residues
        #and a tuple with ('CYS', 'SG', 'disulf') for the disulfide bonds.
        #Separately I create a list composed of tuples containing the resnumbers
        #of the 2 CYS that bound through disulfide bond
        if '_struct_conn.conn_type_id' in cif_dict.keys():
            for i, bound_type in enumerate(cif_dict['_struct_conn.conn_type_id']):

                if bound_type == 'metalc':

                    if cif_dict['_struct_conn.ptnr1_label_comp_id'][i]\
                        in metals\
                        or cif_dict['_struct_conn.ptnr2_label_comp_id'][i]\
                        in metals:

                        if cif_dict['_struct_conn.ptnr1_label_comp_id'][i] in metals:
                            metal_index = '1'
                            res_index = '2'
                        elif cif_dict['_struct_conn.ptnr2_label_comp_id'][i] in metals:
                            metal_index = '2'
                            res_index = '1'

                        substitutions[cif_dict[f'_struct_conn.ptnr{res_index}_auth_seq_id'][i]] =\
                        self.metal_bound_res_parsing(
                        res_auth_seq_id = cif_dict[f'_struct_conn.ptnr{res_index}_auth_seq_id'][i],
                        res_label_atom_id = cif_dict[f'_struct_conn.ptnr{res_index}_label_atom_id'][i],
                        res_label_comp_id = cif_dict[f'_struct_conn.ptnr{res_index}_label_comp_id'][i],
                        metal_label_comp_id = cif_dict[f'_struct_conn.ptnr{metal_index}_label_comp_id'][i]) 
                    
                    else:
                        string = f"Not implemented metal bound, going on pretending nothing happened\n\
                        More info:\n\
                        _struct_conn.ptnr1_label_comp_id = \
                        {cif_dict['_struct_conn.ptnr1_label_comp_id'][i]}\n\
                        _struct_conn.ptnr2_label_comp_id = \
                        {cif_dict['_struct_conn.ptnr2_label_comp_id'][i]}" 

                        print(string)  

                elif bound_type == 'disulf':

                    #if it is present ptnr1_auth_seq_id is the resnumber that remains when transforming it in a pdb
                    try:
                        sulf_bond_tmp = (cif_dict['_struct_conn.ptnr1_auth_seq_id'][i], cif_dict['_struct_conn.ptnr2_auth_seq_id'][i])
                    except:
                        sulf_bond_tmp = (cif_dict['_struct_conn.ptnr1_label_seq_id'][i], cif_dict['_struct_conn.ptnr2_label_seq_id'][i])
                    
                    sulf_bonds.append(sulf_bond_tmp)
                    
                    substitutions[sulf_bond_tmp[0]] = ('CYS', 'SG', 'disulf')
                    substitutions[sulf_bond_tmp[1]] = ('CYS', 'SG', 'disulf')

                else:
                    pass
    
        return substitutions, sulf_bonds, ligand

    def parse_ligands_from_header(self, cif_dict = None, metals = important_lists.metals, trash = important_lists.trash, protein_chain = 'A'):
        """Parses the ligands resnames from a mmcif header"""    

        if type(cif_dict) != dict:
            if  not isinstance(cif_dict, Bio.PDB.MMCIF2Dict.MMCIF2Dict):
                print(not isinstance(cif_dict, Bio.PDB.MMCIF2Dict.MMCIF2Dict))
                raise TypeError(f"Need a dictionary, not a {type(cif_dict)}")
        
        ligand = []

        #record the ligands resname
        for i in cif_dict['_struct_site.details']:
            n = i.split()
            #chooses the right chain only
            if n[-2].strip() == protein_chain:
                n = n[-3]
                n = n.strip()

                if n not in metals:
                    if n not in trash:
                        ligand.append(n)  

        return ligand

    def metal_bound_res_parsing(self,
                                res_auth_seq_id,
                                res_label_atom_id,
                                res_label_comp_id, 
                                metal_label_comp_id):
        """This method is caled by parse_substitutions_PDB method
        returns a tuple with: (resname, binding atom label, metal at which it's bound)"""

        substitution = (res_label_comp_id, res_label_atom_id, metal_label_comp_id)

        return substitution

    def get_ligand_resnum(self, Protein = None, ligand_resnames = None, chain_model_selection = False):
        """Given a Protein instance and a list of Ligand_resnames will return 
        a list containing the ligand resnames and resnumbers
        in order to distinguish ligands with the same resname: [[resname, resnumber], [.., ...], ...]
        
        chain_model_selection :: bool
        if True will only select Protein.model model and Protein.chain chain
        from the given structure"""
        
        if Protein.pdb_file == None or ligand_resnames == None:
            raise TypeError('Need a valid file and ligand_resnames, None is not valid')
        
        # If the ligand resname is a single string I transform it in an iterable object
        ligand_resnames = get_iterable.get_iterable(ligand_resnames)

        if len(ligand_resnames) == 0:
            print("The list of ligands is empty, going on returning a None item")
            return None

        if Protein.file_type == 'cif': 
            p = Bio.PDB.MMCIFParser()

        elif Protein.file_type == 'pdb':
            p = Bio.PDB.PDBParser()

        else:
            raise TypeError(f'"{Protein.file_type}" is not a valid type\n only "cif" and "pdb"')

        struct = p.get_structure(Protein.protein_id, Protein.pdb_file)

        if chain_model_selection == True:
            try:
                #Taking only the right chain and model
                struct = struct[Protein.model][Protein.chain]
            except KeyError:
                struct = struct

        residues = struct.get_residues()

        ligand_list = []

        for residue in residues:
            if residue.resname.strip() in ligand_resnames:
                ligand_list.append([residue.resname.strip(), residue.id[1]])

        return ligand_list

    def get_cysteine_dict(self, Protein = None):
        """Make a dict of the form {resnumber: cysteine_number}
        Where Cysteine number is the number that says if it is the first, second ... Cysteine
        of the seqres

        takes a Protein intance and returns a Protein
        Protein.pdb_file must be a mmCIF file

        Very usefull when you need to know which cysteines make disulf bonds
        but you have to change the resnumbers"""

        if Protein.file_type != "cif":
            raise TypeError(Protein.file_type)
        
        #Iterate on the residues in order to check for cysteines
        cys_dict = {}
        i = 0
        with open(Protein.pdb_file, 'r') as f:
            for line in f:
                if line[0:4] == 'ATOM':

                    if line.split()[5].strip().upper() in important_lists.cyst_resnames:
                        #i is the cysteine_number and line.split()[8].strip() is the residue number
                        if not i in cys_dict.keys():
                            i = i+1
                            cys_dict[line.split()[8].strip()] = str(i)

        Protein.cys_dict = cys_dict

        print(cys_dict)
        print(Protein.sulf_bonds)
        #Checking if everything went the right way
        if Protein.sulf_bonds != None and len(Protein.sulf_bonds) > 0:
            for i in Protein.sulf_bonds:
                for j in i:
                    if j not in Protein.cys_dict.keys():
                        raise ValueError("Something didn't work during the creation of the cys_dict")

        return Protein


    def apply_substitutions(self, protein_id, input_filename, output_filename, substitutions_dict):

        """Applies the substitutions parsed in parse_substitutions_PDB(file_name)
        and returns a Bio.PDB (biopython) structure
            
        substitutions_dict must be a dict (key = residue number, value = new resname)
        """

        import Bio.PDB

        p = Bio.PDB.PDBParser()

        s = p.get_structure(protein_id, input_filename)

        for model in s:
            for chain in model:
                for residue in chain:
                        
                    id = str(residue.id[1])
                    #id = id.replace("'", "").replace("(", "").replace(")", "")
                    id = id.strip()
            
                    if id in substitutions_dict.keys():

                        residue.resname = substitutions_dict[id]

        out = Bio.PDB.PDBIO()
        out.set_structure(s)
        out.save(output_filename)    
        return output_filename


def select_model_chain_custom(Protein = None):
    """Takes a Protein instance containing the filename of a PDB or a mmcif
    Returns a Protein instance with an updated pdb or mmcif file
    using biopython
    selects only a chosen model and chain
    and eliminates disordered atoms chosing only one conformation
    
    protein_chain must be a string
    protein_model must be an integer"""

    import Bio.PDB
    import os

    if Protein == None:
        raise Exception('Protein cannot be None type')

    if Protein.protein_id == None or Protein.pdb_file == None:
        raise ValueError('protein_id and input_filename must be given')
    elif not os.path.exists(Protein.pdb_file):
        raise ValueError(f"{Protein.pdb_file} not found")

    #Bring None arguments to default
    if Protein.model == None:
        Protein.model = 0
    
    if Protein.chain == None:
        Protein.chain = 'A'

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
   
    struct = p.get_structure(Protein.protein_id, Protein.pdb_file)

    if Protein.file_type == 'pdb':
        s = Bio.PDB.PDBIO()

    elif Protein.file_type == 'cif':
        s = Bio.PDB.MMCIFIO()

    s.set_structure(struct[Protein.model][Protein.chain])
    s.save(Protein.pdb_file)

    if not os.path.exists(Protein.pdb_file):
        raise Exception(f'Could not save {Protein.pdb_file} file')

    return Protein


def get_metal_binding_residues_with_no_header(protein_id = None,
                                            pdb_file = None,
                                            mmcif_file = None,
                                            cutoff = 3.0,
                                            substitutions_dict = {},
                                            protein_chain = 'A',
                                            protein_model = 0,
                                            COM_distance = 10.0):

    """This function iterates through the structure many times in order
    to return the metal binding residues through a substitution dictionary

    {residue_id : [residue_name, binding_atom, binding_metal]}

    it can be used both for pdb files and mmcif files (give the path to the files as a string
    to pdb_file or mmcif_file)

    cutoff :: double the maximum distance that a residue's center of mass and a metal ion
    can have to be considered binding default 3.0 angstrom

    if you already have a substitution dictionary and you want to update it give it as input
    as substitutions_dict

    protein_chain :: string default 'A'

    this function is slow and error prone
    and should only be used if there is no mmCIF with a good header
    
    It should not be necessary to change COM_distance because it simply is the distance between the center of mass of
    a residue and the metal that is used to know which atom distances to calculate"""

    if pdb_file == None and mmcif_file == None:
        raise ValueError("I need a pdb_file or a mmcif_file filename cannot both be None type")

    elif pdb_file != None and mmcif_file != None:
        raise ValueError(f"You can only pass a pdb_file or a mmcif_file not both\npdb_file = {pdb_file} mmcif_file = {mmcif_file}")

    elif pdb_file != None:

        if protein_id == None:
            protein_id = pdb_file[0:3] # Get from filename

        p = Bio.PDB.PDBParser()
        structure = p.get_structure(protein_id, pdb_file)

    elif mmcif_file != None:

        if protein_id == None:
            protein_id = mmcif_file[0:3] # Get from filename

        p = Bio.PDB.MMCIFParser()
        structure = p.get_structure(protein_id, mmcif_file)
    
    orient_object = orient.Orient()

    #select the right model
    model = structure[protein_model]
    # select only the right chain
    chain = model[protein_chain.strip().upper()]
    for residue in chain:
        
        if residue.resname.strip().upper() in important_lists.metals:
            for atom in residue:
                    
                #I get a second copy of all the residues in the chain
                #Will have to refactor and clean this mess
                if pdb_file != None:
                    tmp_struct = p.get_structure(protein_id, pdb_file)
                else:
                    tmp_struct = p.get_structure(protein_id, mmcif_file)
                tmp_struct = tmp_struct[protein_model][protein_chain.strip().upper()]
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
                                substitutions_dict[str(other_residue.id[1])] = [other_residue.resname.strip().upper(), TMP_atom_dist[1], residue.resname.strip().upper()]

        
    #useless variables
    COM_1 = None
    COM_2 = None

    return substitutions_dict


def get_disulf_bonds_with_no_header(protein_id = None,
                                    pdb_file = None,
                                    mmcif_file = None,
                                    cutoff = 3.0,
                                    substitutions_dict = {},
                                    sulf_bonds = [],
                                    protein_chain = 'A',
                                    protein_model = 0):
    """This function iterates through the structure many times in order
    to return the disulf bonds through a substitution dictionary and a list of the binded couples

    {residue_id : [residue_name, binding_atom, binding_metal]}
    and
    [(cys_id, cys_id), (cys_id, ...), ...]

    return substitutions_dict, sulf_bonds

    it can be used both for pdb files and mmcif files (give the path to the files as a string
    to pdb_file or mmcif_file)

    cutoff :: double the maximum distance that two CYS S atoms
    can have to be considered binding default 3.0 angstrom

    if you already have a substitution dictionary and you want to update it give it as input
    as substitutions_dict

    protein_chain :: string default 'A'

    this function is slow and error prone
    and should only be used if there is no mmCIF with a good header"""

    if pdb_file == None and mmcif_file == None:
        raise ValueError("I need a pdb_file or a mmcif_file filename cannot both be None type")

    elif pdb_file != None and mmcif_file != None:
        raise ValueError(f"You can only pass a pdb_file or a mmcif_file not both\npdb_file = {pdb_file} mmcif_file = {mmcif_file}")

    elif pdb_file != None:

        if protein_id == None:
            protein_id = pdb_file[0:3] # Get from filename

        p = Bio.PDB.PDBParser()
        structure = p.get_structure(protein_id, pdb_file)

    elif mmcif_file != None:

        if protein_id == None:
            protein_id = mmcif_file[0:3] # Get from filename

        p = Bio.PDB.MMCIFParser()
        structure = p.get_structure(protein_id, mmcif_file)

    #select the right model
    model = structure[protein_model]
    # select only the right chain
    chain = model[protein_chain.strip().upper()]
    for residue in chain:
        
        if residue.resname.strip().upper() in important_lists.cyst_resnames:

            #I get a second copy of all the residues in the chain
            all_residues = chain.get_residues()
            #and iterate though them
            for other_residue in all_residues:
                
                if other_residue.resname.strip().upper() in important_lists.cyst_resnames:
                    #I avoid scanning CYS against its self and agaisnt the already scanned ones
                    if other_residue.id[1] > residue.id[1]:

                        distance = (residue['SG'].coord[0] - other_residue['SG'].coord[0])**2. + (residue['SG'].coord[1] - other_residue['SG'].coord[1])**2. + (residue['SG'].coord[2] - other_residue['SG'].coord[2])**2.
                        distance = distance ** (0.5)

                        if distance <= cutoff:

                            substitutions_dict[str(residue.id[1])] = ['CYS', 'SG', 'disulf']
                            substitutions_dict[str(other_residue.id[1])] = ['CYS', 'SG', 'disulf']

                            sulf_bonds.append((str(residue.id[1]), str(other_residue.id[1])))

    return substitutions_dict, sulf_bonds

def get_organic_ligands_with_no_header(protein_id = None,
                                    pdb_file = None,
                                    mmcif_file = None,
                                    protein_chain = 'A',
                                    protein_model = 0):
    """This function iterates through the structure to get the organic ligand
    returning a list of lists containing
    [[resname, resnumber], [resname, resnumber], ...]

    it can be used both for pdb files and mmcif files (give the path to the files as a string
    to pdb_file or mmcif_file)

    protein_chain :: string default 'A'

    this function is slow and error prone
    and should only be used if there is no mmCIF with a good header"""

    if pdb_file == None and mmcif_file == None:
        raise ValueError("I need a pdb_file or a mmcif_file filename cannot both be None type")

    elif pdb_file != None and mmcif_file != None:
        raise ValueError(f"You can only pass a pdb_file or a mmcif_file not both\npdb_file = {pdb_file} mmcif_file = {mmcif_file}")

    elif pdb_file != None:

        if protein_id == None:
            protein_id = pdb_file[0:3] # Get from filename

        p = Bio.PDB.PDBParser()
        structure = p.get_structure(protein_id, pdb_file)

    elif mmcif_file != None:

        if protein_id == None:
            protein_id = mmcif_file[0:3] # Get from filename

        p = Bio.PDB.MMCIFParser()
        structure = p.get_structure(protein_id, mmcif_file)

    ligand_list = []

    #select the right model
    model = structure[protein_model]
    # select only the right chain
    chain = model[protein_chain.strip().upper()]

    for residue in chain:

        #it's an hetero atom
        if residue.id[0].strip() != '':

            if residue.resname.strip().upper() not in important_lists.metals and residue.resname.strip().upper() not in important_lists.trash:

                ligand_list.append([residue.resname.strip().upper(), str(residue.id[1])])

    return ligand_list

def remove_trash_metal_ions(pdb_file):
    """This function removes unwanted metal ions
    that are still inside the structure after it went through prody
    selection

    This is a brutal function I will need to do a better job
    the file must be a pdb"""

    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    with open(pdb_file, 'w') as w:
        for line in lines:

            if line[0:4] == 'ATOM' or line[0:6] == 'HETATM' or line[0:3] == 'TER':

                if line[17:20].strip().upper() not in important_lists.trash_ions:

                    w.write(f"{line.strip()}\n")

    return pdb_file