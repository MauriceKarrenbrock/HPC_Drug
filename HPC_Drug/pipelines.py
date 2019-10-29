# This file contains the possible pipelines depending from input

import file_manipulation
import structures

def choose_pipeline(*args, **kwargs):
    """Takes the parsed program input
    and returns the right pipeline class"""
    
    input_dict = kwargs

    if input_dict['protein_filetype'] == 'pdb':
        # Protein is a pdb file
        if (input_dict['ligand'] == None or input_dict['ligand'] == '$protein'):
            #ligand shall be taken from pdb (at your own risk)
            
            return PdbNoLigand_Pipeline(protein = input_dict['protein'],\
                protein_filetype = input_dict['protein_filetype'],\
                    local = input_dict['local'], filepath = input_dict['filepath'], ligand = input_dict['ligand'],\
                        ligand_elaboration_program = input_dict['ligand_elaboration_program'],\
                            PDB_model = input_dict['PDB_model'])

        else:
            # ligand is given

            return PdbLigand_Pipeline()

    elif input_dict['protein_filetype'] == 'cif':
        # Protein is a pdb file
        if (input_dict['ligand'] == None or input_dict['ligand'] == '$protein'):
            #ligand shall be taken from pdb (at your own risk)
            
            return CifNoLigand_Pipeline()

        else:
            # ligand is given

            return CifLigand_Pipeline()
    else:
        raise NotImplementedError('The pipeline for {} file type was not implemented'.format(input_dict['protein_filetype']))

class Pipeline(object):
    """General pipeline class"""
    
    def __init__(self, protein = None, protein_filetype = None,\
        local = None, filepath = None, ligand = None, ligand_elaboration_program = None, PDB_model = None):
        
        self.protein_id = protein
        self.protein_filetype = protein_filetype
        self.local = local
        self.protein_filename = filepath
        self.model = PDB_model
        
        #ligand can be a pdb file or a smiles
        self.ligand_filename = ligand
        self.ligand_elaboration_program = ligand_elaboration_program

    def download(self):
        return file_manipulation.download_protein_structure(self.protein_id, self.protein_filetype, None)


class CifNoLigand_Pipeline(Pipeline):
    """Protein is given as a PDBx/mmCIF
    The ligand is not given as input
    Will be searched in the protein file"""
    pass

class CifLigand_Pipeline(Pipeline):
    """Protein is given as a PDBx/mmCIF
    The ligand is given as input"""
    pass


class PdbNoLigand_Pipeline(Pipeline):
    """Protein is given as a PDB
    The ligand is not given as input
    Will be searched in the protein file"""
    
    def execute(self):
        """The execution of the pipeline"""
        import os

        #If requested in input will download pdb file
        #If the given local file doesn't exist will download pdb
        if self.local == 'no' or self.local == None or (not os.path.exists(self.protein_filename)):
            self.protein_filename = self.download()
        elif self.local != 'yes':
            raise ValueError(f'local cannot be {self.local} it can only be "yes" "no" or omitted')
        else:
            #protein filename points to an existing file
            #nothing is done
            pass

        cruncer = file_manipulation.ProteinCruncer(self.protein_filetype)
        Protein = cruncer.get_protein(self.protein_id, self.protein_filename, self.model, None)
        Ligand = cruncer.get_ligand(self.protein_id, self.protein_filename, self.model, None)

        Protein.write_PDB(f"{self.protein_id}_protein.pdb")
        Ligand.write_PDB(f"{self.protein_id}_ligand.pdb")





class PdbLigand_Pipeline(Pipeline):
    """Protein is given as a PDB
    The ligand is given as input"""
    pass