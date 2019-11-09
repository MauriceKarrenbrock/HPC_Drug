# This file contains the classes and functions to create structures
import file_manipulation
class Structure(object):
    """Raises an Exception for not implemented structure type"""
    def __init__(self):
        raise NotImplementedError('This structure was not implemented yet')
        
class Protein(Structure):
    """The general Protein class"""

    def __init__(self, protein_id = None, filename = None,
                structure = None, substitutions_dict = None,
                sulf_bonds = None, seqres = None,
                file_type = 'cif', model = None):
        self.protein_id = protein_id
        self.model = model
        self.file_type = file_type
        self.filename = filename

        if protein_id == None:
            raise Exception('need a protein id')
        elif self.filename == None:
            self.filename = protein_id + '.' + file_type
        
        self.structure = structure

        self.substitutions_dict = substitutions_dict
        self.sulf_bonds = sulf_bonds
        self.seqres = seqres
    
    def write_PDB(self, filename = None):
        import file_manipulation as fm

        if filename == None:
            filename = self.protein_id + '_protein.pbd'

        w = fm.ProteinCruncer('pdb')
        w.write_PDB(self.structure, filename)


class Ligand(Structure):
    """The general Ligand class"""

    def __init__(self, ligand_resname = None, filename = None, structure = None, file_type = 'pdb'):
        self.ligand_resname = ligand_resname
        self.file_type = file_type

        self.filename = filename

        self.structure = structure
    
    def write_PDB(self, filename):
        import file_manipulation as fm

        if filename == None:
            filename = self.ligand_resname + '_ligand.pbd'

        w = fm.ProteinCruncer('pdb')
        w.write_PDB(self.structure, filename)

        return filename
