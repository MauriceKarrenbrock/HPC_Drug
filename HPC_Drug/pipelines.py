# This file contains the possible pipelines depending from input

import subprocess
import file_manipulation
import structures
import pipeline_functions
import funcs4primadorac

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
                                    ligand_elaboration_program_path = input_dict['ligand_elaboration_program_path'],
                                    Protein_model = input_dict['Protein_model'],
                                    Protein_chain = input_dict['Protein_chain'],
                                    ph = input_dict['ph'],
                                    repairing_method = input_dict['repairing_method'],
                                    MD_program = input_dict['MD_program'],
                                    MD_program_path = input_dict['MD_program_path'])

        else:
            # ligand is given

            return Ligand_Pipeline()
        
    else:
        raise NotImplementedError('The pipeline for {} file type was not implemented'.format(input_dict['protein_filetype']))

class Pipeline(object):
    """General pipeline class"""
    
    def __init__(self,
                protein = None,
                protein_filetype = None,
                local = None,
                filepath = None,
                ligand = None,
                ligand_elaboration_program = 'primadorac',
                ligand_elaboration_program_path = 'primadorac.bash',
                Protein_model = None,
                Protein_chain = None,
                ph = 7.0,
                repairing_method = 'pdbfixer',
                MD_program = None,
                MD_program_path = None):
        
        self.protein_id = protein
        self.protein_filetype = protein_filetype
        self.local = local
        self.protein_filename = filepath
        self.model = Protein_model
        self.chain = Protein_chain

        self.ph = ph
        self.repairing_method = repairing_method
        
        #ligand can be a pdb file or a smiles
        self.ligand_filename = ligand
        self.ligand_elaboration_program = ligand_elaboration_program
        self.ligand_elaboration_program_path = ligand_elaboration_program_path

        self.MD_program = MD_program
        self.MD_program_path = MD_program_path

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
        

        #Declaring protein instance
        Protein = structures.Protein(protein_id = self.protein_id,
                                    filename = self.protein_filename,
                                    model = self.model, chain = self.chain)
        


        #Parses the mmcif for sulf bonds and organic ligands
        #takes information about the residues binding metals
        #and records seqres
        #ligand_resnames is a list containing the ligands resnames
        Protein, ligand_resnames = pipeline_functions.parse(Protein)
        
        repairer = file_manipulation.PDBRepair()
        
        #adds missing atoms and residues
        #changes non standard residues to standard ones
        Protein.protein_filename = repairer.add_missing_atoms(Protein.protein_id,
                                                            Protein.filename, self.repairing_method,
                                                            None, ph = self.ph, add_H = False)

        #Protein.protein_filename = subst_parser.apply_substitutions(self.protein_id, self.protein_filename, self.protein_filename, substitutions_dict)

        cruncer = file_manipulation.ProteinCruncer(Protein.file_type)
        
        #gets the protein's structure from the pdb
        #The only HETATM remaining are the metal ions
        Protein.structure = cruncer.get_protein(protein_id = Protein.protein_id,
                                                filename =Protein.filename,
                                                structure = None)

        Protein.filename = Protein.write_PDB(f"{Protein.protein_id}_protein.pdb")
        
        #extracts the ligands structures (if any) from the pdb
        #creates a ligand structure and writes a pdb
        #for any found organic ligand
        if len(ligand_resnames) != 0:
            #Ligand from now on is a list of Ligand objects
            Ligand = []
            for i, resname in enumerate(pipeline_functions.get_iterable(ligand_resnames)):
                Ligand.append(None)

                Ligand[i] = structures.Ligand(resname, None, None, 'pdb')
                Ligand[i].structure = cruncer.get_ligand(Protein.protein_id,
                                                        Protein.protein_filename,
                                                        resname,
                                                        None)

                print(Ligand[i].structure)

                #tries to write the ligand pdb
                try:
                    Ligand[i].ligand_filename = Ligand[i].write_PDB(f'{Protein.protein_id}_{resname}_ligand{i}.pdb')
                except TypeError as err:
                    raise TypeError(f'{err.args}\ncannot make ligand pdb file for {resname}')
                except Exception as err:
                    raise Exception(f'{err} {err.args}\ncannot make ligand pdb file for {resname}')
        else:
            Ligand = None

        #use primadorac to get topology and parameter files for any given ligand
        if self.ligand_elaboration_program == 'primadorac':
        
            #runs primadorac and returns the Ligand list updated with the prm, tpg and new pdb files
            Ligand = funcs4primadorac.run_primadorac(ligand_list = Ligand,
                                                    primadorac_path = self.ligand_elaboration_program_path,
                                                    ph = self.ph)


        else:
            raise NotImplementedError(self.ligand_elaboration_program)

        if self.ligand_elaboration_program == None or self.ligand_elaboration_program_path == None:
            raise Exception('Need a ligand elaboration program and a path')
        
        elif self.ligand_elaboration_program == 'primadorac':
            #DA FARE
            pass
        
        else:
            raise NotImplementedError(self.ligand_elaboration_program)

        if self.MD_program == None or self.MD_program_path == None:
            raise Exception('Need a MD program and a path')

        elif self.MD_program == 'gromacs':
            import funcs4gromacs

            Protein = funcs4gromacs.residue_substitution(Protein, 'standard')
            Protein = pipeline_functions.get_seqres_PDB(Protein)

        elif self.MD_program == 'orac':
            import funcs4orac

            Protein = funcs4orac.residue_substitution(Protein, 'standard')
            Protein = pipeline_functions.get_seqres_PDB(Protein)

        else:
            raise NotImplementedError(self.MD_program)
        
        


class Ligand_Pipeline(Pipeline):
    """Protein is given as a mmcif
    The ligand is given as input"""
    pass
