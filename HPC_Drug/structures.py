# This file contains the classes and functions to create structures
from HPC_Drug import file_manipulation
class Structure(object):
    """Raises an Exception for not implemented structure type"""
    def __init__(self):
        raise NotImplementedError('This structure was not implemented yet')
        
class Protein(Structure):
    """The general Protein class"""

    def __init__(self,
                protein_id = None,
                filename = None,
                structure = None,
                substitutions_dict = None,
                sulf_bonds = None,
                seqres = None,
                file_type = 'cif',
                model = None,
                chain = None,
                cys_dict = None,
                gro_file = None,
                top_file = None):
                
        self.protein_id = protein_id
        self.model = model
        self.chain = chain
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

        self.cys_dict = cys_dict

        self.gro_file = gro_file
        self.top_file = top_file
    
    def write_PDB(self, filename = None, struct_type = 'prody'):
        from HPC_Drug import file_manipulation as fm

        if filename == None:
            filename = self.protein_id + '_protein.pbd'

        if struct_type == 'prody':

            w = fm.ProteinCruncer('pdb')
            w.write_PDB(self.structure, filename)

        elif struct_type == 'biopython':
            import Bio.PDB 

            p = Bio.PDB.PDBIO()
            p.set_structure(self.structure)
            p.save(filename)

        return filename


class Ligand(Structure):
    """The general Ligand class"""

    def __init__(self,
                ligand_resname = None,
                filename = None,
                structure = None,
                file_type = 'pdb',
                topology_file = None,
                param_file = None,
                res_number = None,
                itp_file = None):
        self.ligand_resname = ligand_resname
        self.file_type = file_type

        self.filename = filename

        self.structure = structure

        self.topology_file = topology_file
        self.param_file = param_file

        self.res_number = res_number

        self.itp_file = itp_file
    
    def write_PDB(self, filename = None, struct_type = 'prody'):
        from HPC_Drug import file_manipulation as fm

        if filename == None:
            filename = self.ligand_resname + '_lgand.pbd'

        
        if struct_type == 'prody':

            w = fm.ProteinCruncer('pdb')
            w.write_PDB(self.structure, filename)

        elif struct_type == 'biopython':
            import Bio.PDB 

            p = Bio.PDB.PDBIO()
            p.set_structure(self.structure)
            p.save(filename)


        return filename
