# This file contains the possible pipelines depending from input

import file_manipulation
import structures

def choose_pipeline(*args, **kwargs):
    """Takes the parsed program input
    and returns the right pipeline class"""
    
    input_dict = kwargs

    if input_dict['protein_filetype'] == 'cif':
        # Protein is a pdb file
        if (input_dict['ligand_in_protein'] == None or input_dict['ligand_in_protein'] == 'yes'):
            #ligand will be taken from protein file (at your own risk)
            
            return NoLigand_Pipeline(protein = input_dict['protein'],
                                    protein_filetype = input_dict['protein_filetype'],
                                    local = input_dict['local'],
                                    filepath = input_dict['filepath'],
                                    ligand = input_dict['ligand'],
                                    ligand_elaboration_program = input_dict['ligand_elaboration_program'],
                                    Protein_model = input_dict['Protein_model'],
                                    ph = input_dict['ph'],
                                    repairing_method = input_dict['repairing_method'])

        else:
            # ligand is given

            return Ligand_Pipeline()
        
    else:
        raise NotImplementedError('The pipeline for {} file type was not implemented'.format(input_dict['protein_filetype']))

class Pipeline(object):
    """General pipeline class"""
    
    def __init__(self, protein = None, protein_filetype = None,
                local = None, filepath = None, ligand = None,
                ligand_elaboration_program = None,
                Protein_model = None, ph = 7.0,
                repairing_method = 'pdbfixer'):
        
        self.protein_id = protein
        self.protein_filetype = protein_filetype
        self.local = local
        self.protein_filename = filepath
        self.model = Protein_model

        self.ph = ph
        self.repairing_method = repairing_method
        
        #ligand can be a pdb file or a smiles
        self.ligand_filename = ligand
        self.ligand_elaboration_program = ligand_elaboration_program

    def download(self):
        return file_manipulation.download_protein_structure(self.protein_id, self.protein_filetype, None)


class NoLigand_Pipeline(Pipeline):
    """Protein is given as a mmcif
    The ligand is not given as input
    Will be searched in the protein file"""
    
    def execute(self):
        """The execution of the pipeline"""
        import os

        #If requested in input will download pdb file
        #If the given local file doesn't exist will download pdb
        if self.local == 'no' or self.local == None or (not os.path.exists(self.protein_filename)):
            self.protein_filename = self.download()
            self.protein_pdb = file_manipulation.download_protein_structure(self.protein_id,
                                                                            self.protein_filetype,
                                                                            None)
        elif self.local != 'yes':
            raise ValueError(f'local cannot be {self.local} it can only be "yes" "no" or omitted')
        else:
            #protein filename points to an existing file
            #nothing is done
            pass
        

        #Declaring protein and ligand instances
        Protein = structures.Protein(protein_id = self.protein_id,
                                    filename = self.protein_filename,
                                    model = self.model)
        
        Ligand = structures.Ligand()

        #Parsing substitutions, sulf bonds and the resname of the organic ligand
        subst_parser = file_manipulation.SubstitutionParser()

        Protein.substitutions_dict, Protein.sulf_bonds, Ligand.ligand_resname = subst_parser.parse_substitutions_PDB(self.protein_filename)
        
        #converting the mmcif file in a pdb
        #because pdbfixer has some issues with mmcif

        Protein.filename = file_manipulation.mmcif2pdb_custom(Protein.protein_id, Protein.filename, Protein.model)
        Protein.file_type = 'pdb'

        
        repairer = file_manipulation.PDBRepair()
        
        #adds missing atoms and residues
        #changes non standard residues to standard ones
        #adds hydrogens according to given ph (default 7.0)
        Protein.protein_filename = repairer.add_missing_atoms(Protein.protein_id,
                                                            Protein.filename, self.repairing_method,
                                                            None, ph = self.ph)

        #Protein.protein_filename = subst_parser.apply_substitutions(self.protein_id, self.protein_filename, self.protein_filename, substitutions_dict)

        cruncer = file_manipulation.ProteinCruncer(Protein.file_type)
        
        Protein.structure = cruncer.get_protein(Protein.protein_id,
                                    Protein.filename, Protein.model, None)

        Protein.write_PDB(f"{Protein.protein_id}_protein.pdb")
        
        if len(Ligand.ligand_resname) != 0:
            for ligand in Ligand.ligand_resname:

                Ligand.structure = cruncer.get_ligand(Protein.protein_id,
                                            Protein.protein_filename, ligand,
                                            Protein.model, None)

                print(Ligand.structure)

                #tries to write the ligand pdb
                try:
                    temp = Ligand.write_PDB(Protein.protein_id + '_' + ligand +'_ligand.pdb')
                except TypeError as err:
                    raise TypeError(f'{err.args}\ncannot make ligand pdb file for {ligand}')
                except Exception as err:
                    raise Exception(f'{err} {err.args}\ncannot make ligand pdb file for {ligand}')


class Ligand_Pipeline(Pipeline):
    """Protein is given as a mmcif
    The ligand is given as input"""
    pass

class Preprocessing_Pipeline(Pipeline):
    """Contains the methods to clean, parse and split the
    mmcif file in a pdb file for the protein and the metallic
    atoms, and another for each organic ligand (if found)"""

    def __init__(self):
        pass

    def parse(self, protein):
        """Takes a Protein object that must have a valid Protein.filename
        Parses the protein file
        for sostitutions, conversts the mmCIF to a PDB file
        and returns a Protein object"""
        pass