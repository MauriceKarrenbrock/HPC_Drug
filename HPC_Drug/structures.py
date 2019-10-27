# This file contains the classes and functions to create structures
class Structure(object):
    """Raises an Exception for not implemented structure type"""
    def __init__(self):
        raise Exception(NotImplementedError)
        
class Protein(Structure):
    """Gets the structure of one module from a PDBx/mmCIF or PDB file"""

    import Bio.PDB

    def __init__(self, protein_id, filename = None, file_type = 'cif', model = 0):
        self.structure_file = filename
        self.protein_id = protein_id
        self.model = model
        self.file_type = file_type
        if filename == None:
            self.structure_file = protein_id + '.cif'
        self.structure = self.create_model_structure(self.protein_id, self.structure_file, self.model)

    def download_protein_structure(protein_id, file_type = 'mmCif', pdir = None):
        import Bio.PDB
        print(file_type)

        if file_type == 'cif':
            file_type == 'mmCif'
        file_type = None
        print(file_type)

        pdbl = Bio.PDB.PDBList()
        pdbl.retrieve_pdb_file(protein_id, file_type, pdir)
    
    def create_structure(self, protein_id, filename):
        """ create_structure(protein_id, filename)
        First parses the .cif file
        then gets the  general protein structure
        then gets module number "model" (default model = 0) to be sure to have only one module
        returns the module "model" structure as protein stucture"""

        if self.file_type == 'cif':
            import Bio.PDB.MMCIFParser
            parser = Bio.PDB.MMCIFParser()
        else:
            parser = Bio.PDB.PDBParser()
        
        structure = parser.get_structure(protein_id, filename)
        
        return structure

    def create_model_structure(self, protein_id, filename, model = 0):
        
        model_structure = self.create_structure(protein_id, filename)
        model_structure = model_structure[model]

        return model_structure

    def create_ligand_structure(self, protein_id, filename):
        pass

class Ligand(Structure):
    pass

class LigandFromProteinFile(Ligand):
    pass

class LigandFromInput(Ligand):
    #mi limito a registrare l'input poi do tutto a primadorac
    pass