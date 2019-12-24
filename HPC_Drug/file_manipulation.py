# Contains the classes for file manipulation like PDB & mmCIF

from HPC_Drug import pipeline_functions
from HPC_Drug import structures
from HPC_Drug import important_lists
import Bio.PDB
import Bio.PDB.MMCIF2Dict
import os
import prody





def download_protein_structure(protein_id, file_type = 'cif', pdir = None):
    """The function downloads a PDB or a mmCIF from wwwPDB in a selected directory
    the default directory is the working directory
    it returns the filename (str)"""
    

    if file_type == 'cif':
        file_type = 'mmCif'
    
    elif file_type != 'pdb' and file_type != 'cif':
        raise ValueError(f"Must be 'pdb' or 'cif' not {file_type}")
    
    if pdir == None:
        pdir = os.getcwd()

    pdbl = Bio.PDB.PDBList()
    filename = pdbl.retrieve_pdb_file(protein_id, False, pdir, file_format = file_type, overwrite = True)
    
    if not os.path.exists(filename):
        raise FileNotFoundError(f'Was not able to download the protein or to find {filename}')
    
    else:
        return filename


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
                        repairing_method = 'pdbfixer',
                        output_filename = None, ph = 7.0, add_H = False):
        """This function returns a repaired PDBx/mmCIF with missing atoms but not hydrogens
        chooses the right repairing method with repairing_method (pdbfixer is default)
        The returned PDB is called output_filename
        if output_filename is None (default) the file is saved as
        {protein_id}_repaired.cif in the working directory
        returns the output_filename

        does only work on the protein, not on the ligand
        
        The input file MUST be a PDBx/mmCIF"""

        if repairing_method == 'pdbfixer':
            import pdbfixer
            import simtk.openmm.app
            import os

            fixer = self.fix_pdbfixer(input_filename, ph = ph)

            if output_filename == None:
                output_filename = pdb_id + '_repaired.cif'
            
            with open(output_filename, 'w') as f:
                simtk.openmm.app.pdbxfile.PDBxFile.writeFile(fixer.topology, fixer.positions, f, keepIds= True)
            
            output_filename = os.getcwd() + '/' + output_filename
            if not os.path.exists(output_filename):
                raise Exception(f'was not able to write {output_filename}')
            else:
                
                
                return output_filename
        
        else:
            raise NotImplementedError(f'{repairing_method}')


    def fix_pdbfixer(self, input_filename, ph = 7.0, add_H = False):
        """This method is called by add_missing_atoms
        Reads an PDBx/mmCIF"""
        
        import pdbfixer
        import simtk.openmm.app

        with open(input_filename, 'r') as f:
            fixer = pdbfixer.PDBFixer(pdbxfile = f)
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        #fixer.removeHeterogens(False)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        #fixer.addMissingHydrogens(8.0)
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
            cif_dict[key] = pipeline_functions.get_iterable(cif_dict[key])

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
        
        if Protein.filename == None or ligand_resnames == None:
            raise TypeError('Need a valid mmcif and ligand_resnames, None is not valid')

        if len(ligand_resnames) == 0:
            print("The list of ligands is empty, going on returning a None item")
            return None

        if Protein.file_type == 'cif': 
            p = Bio.PDB.MMCIFParser()

        elif Protein.file_type == 'pdb':
            p = Bio.PDB.PDBParser()

        else:
            raise TypeError(f'"{Protein.file_type}" is not a valid type\n only "cif" and "pdb"')

        struct = p.get_structure(Protein.protein_id, Protein.filename)

        if chain_model_selection == True:
            #Taking only the right chain and model
            struct = struct[Protein.model][Protein.chain]

        residues = struct.get_residues()

        ligand_list = []

        for residue in residues:
            if residue.resname in ligand_resnames:
                ligand_list.append([residue.resname, residue.id[1]])

        return ligand_list

    def get_cysteine_dict(self, Protein = None):
        """Make a dict of the form {resnumber: cysteine_number}
        Where Cysteine number is the number that says if it is the first, second ... Cysteine
        of the seqres

        takes a Protein intance and returns a Protein
        Protein.filename must be a mmCIF file

        Very usefull whern you need to know which cysteines make disulf bonds
        but you have to change the resnumbers"""

        if Protein.file_type != "cif":
            raise TypeError(Protein.file_type)
        
        #Iterate on the residues in order to check for cysteines
        cys_dict = {}
        i = 0
        with open(Protein.filename, 'r') as f:
            for line in f:
                if "ATOM" in line:
                    if line[0] == "A":

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
    Returns a Protein instance with an updated pdb file
    using biopython
    selects only a chosen model and chain
    and eliminates disordered atoms chosing only one conformation
    
    protein_chain must be a string
    protein_model must be an integer"""

    import Bio.PDB
    import os

    if Protein == None:
        raise Exception('Protein cannot be None type')

    if Protein.protein_id == None or Protein.filename == None:
        raise ValueError('protein_id and input_filename must be given')
    elif not os.path.exists(Protein.filename):
        raise ValueError(f"{Protein.filename} not found")

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
   
    struct = p.get_structure(Protein.protein_id, Protein.filename)

    s = Bio.PDB.PDBIO()
    s.set_structure(struct[Protein.model][Protein.chain])
    s.save(Protein.filename)

    if not os.path.exists(Protein.filename):
        raise Exception(f'Could not save {Protein.filename} file')

    return Protein
