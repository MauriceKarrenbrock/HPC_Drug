class Structure(object):
    """Raises an Exception for not implemented structure type"""
    def __init__(self):
        raise Exception(NotImplementedError)

class Protein_Cif(Structure):
    """Gets the structure of one module from a PDBx/mmCIF file"""

    def __init__(self, protein_id, filename = None, model = 0):
        if filename != None:
            self.structure_file = filename
        else:
            self.structure_file = protein_id + '.cif'
        
        self.structure = self.create_structure(protein_id, self.structure_file, model)
    
    def create_structure(self, protein_id, filename, model = 0):
        """ create_structure(protein_id, filename)
        First parses the .cif file
        then gets the  general protein structure
        then gets module number "model" (default model = 0) to be sure to have only one module
        returns the module "model" structure as protein stucture"""

        import Bio.PDB
        import Bio.PDB.MMCIFParser
        
        parser = Bio.PDB.MMCIFParser()
        structure = parser.get_structure(protein_id, filename)
        model_structure = structure[model]

        return model_structure
