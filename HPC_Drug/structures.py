# This file contains the classes and functions to create structures
import file_manipulation
class Structure(object):
    """Raises an Exception for not implemented structure type"""
    def __init__(self):
        raise NotImplementedError('This structure was not implemented yet')
        
class Protein(Structure):
    """The general Protein class"""

    def __init__(self, protein_id = None, filename = None, structure = None, file_type = None, model = None):
        self.protein_id = protein_id
        self.model = model
        self.file_type = file_type
        if filename == None:
            self.structure_file = protein_id + '.' + file_type
        self.structure = structure
    
    def write_PDB(self, filename = None):
        import file_manipulation as fm

        if filename == None:
            filename = self.protein_id + '_protein.pbd'

        w = fm.ProteinCruncer('pdb')
        w.write_PDB(self.structure, filename)


class Ligand(Structure):
    """The general Ligand class"""

    def __init__(self, ligand_id = None, filename = None, structure = None, file_type = None):
        self.ligand_id = ligand_id
        self.file_type = file_type
        if filename == None:
            self.structure_file = ligand_id + '.' + file_type
        self.structure = structure
    
    def write_PDB(self, filename = None):
        import file_manipulation as fm

        if filename == None:
            filename = self.ligand_id + '_ligand.pbd'

        w = fm.ProteinCruncer('pdb')
        w.write_PDB(self.structure, filename)
