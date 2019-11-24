# This file contains the possible pipelines depending from input

import subprocess
from HPC_Drug import file_manipulation
from HPC_Drug import structures
from HPC_Drug import pipeline_functions
from HPC_Drug import funcs4primadorac
from HPC_Drug import funcs4gromacs
from HPC_Drug import funcs4orac

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
                                    MD_program_path = input_dict['MD_program_path'],
                                    protein_prm_file = input_dict['protein_prm_file'],
                                    protein_tpg_file = input_dict['protein_tpg_file'],
                                    solvent_pdb = input_dict['solvent_pdb'])

        else:
            # ligand is given

            return Ligand_Pipeline()
        
    else:
        raise NotImplementedError('The pipeline for {} file type was not implemented'.format(input_dict['protein_filetype']))

class Pipeline(object):
    """General pipeline class"""
    
    def __init__(self,
                protein = None,
                protein_filetype = 'cif',
                local = None,
                filepath = None,
                ligand = None,
                ligand_elaboration_program = 'primadorac',
                ligand_elaboration_program_path = 'primadorac.bash',
                Protein_model = 0,
                Protein_chain = 'A',
                ph = 7.0,
                repairing_method = 'pdbfixer',
                MD_program = None,
                MD_program_path = None,
                protein_prm_file = None,
                protein_tpg_file = None,
                solvent_pdb = None):
        
        self.protein_id = protein
        self.protein_filetype = protein_filetype

        if self.protein_filetype == None:
            self.protein_filetype = 'cif'

        self.local = local
        self.protein_filename = filepath

        self.model = Protein_model
        if self.model == None:
            self.model = 0
        elif type(self.model) == str:
            self.model = int(self.model)

        self.chain = Protein_chain
        if self.chain == None:
            self.chain = 'A'
        else:
            self.chain = self.chain.upper()

        self.ph = ph
        if type(self.ph) == str:
            self.ph = float(self.ph)
         
        self.repairing_method = repairing_method
        
        #ligand can be a pdb file or a smiles
        self.ligand_filename = ligand
        self.ligand_elaboration_program = ligand_elaboration_program
        self.ligand_elaboration_program_path = ligand_elaboration_program_path

        self.MD_program = MD_program
        self.MD_program_path = MD_program_path

        self.protein_prm_file = protein_prm_file
        self.protein_tpg_file = protein_tpg_file
        self.solvent_pdb = solvent_pdb

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
        #ligand_resnames is a list containing the ligands
        # resnames and resnumbers [[resname, resnumber], [..., ...], ...]
        Protein, ligand_resnames = pipeline_functions.parse(Protein)
        
        repairer = file_manipulation.PDBRepair()
        
        #adds missing atoms and residues
        #changes non standard residues to standard ones
        #takes a PDBx/mmCIF and returns a PDB
        Protein.filename = repairer.add_missing_atoms(Protein.protein_id,
                                                            Protein.filename, self.repairing_method,
                                                            None, ph = self.ph, add_H = False)

        Protein.file_type = 'pdb'

        #selects only a selected model and chain, and keeps only one conformation for any disordered atom
        Protein = file_manipulation.select_model_chain_custom(Protein = Protein)

        cruncer = file_manipulation.ProteinCruncer(Protein.file_type)
        
        #gets the protein's structure from the pdb
        #The only HETATM remaining are the metal ions
        Protein.structure = cruncer.get_protein(protein_id = Protein.protein_id,
                                                filename = Protein.filename,
                                                structure = None)

        #Temporarely save th filename that contains the ligand in order to extract it later
        tmp_pdb_file = Protein.filename

        #Write Protein only pdb
        Protein.filename = Protein.write_PDB(f"{Protein.protein_id}_protein.pdb")
        
        #extracts the ligands structures (if any) from the pdb
        #creates a ligand structure and writes a pdb
        #for any found organic ligand
        if ligand_resnames == None:
            Ligand = None
        
        elif len(ligand_resnames) != 0:
            #Ligand from now on is a list of Ligand objects
            Ligand = []

            #ligand res is: [resname, resnumber]
            for i, ligand_res in enumerate(pipeline_functions.get_iterable(ligand_resnames)):
                Ligand.append(None)

                Ligand[i] = structures.Ligand(ligand_resname = ligand_res[0],
                                            filename = None,
                                            structure = None,
                                            file_type = 'pdb',
                                            topology_file = None,
                                            param_file = None,
                                            res_number= ligand_res[1])
                                            
                Ligand[i].structure = cruncer.get_ligand(protein_id = Protein.protein_id,
                                                        filename = tmp_pdb_file,
                                                        ligand_name = Ligand[i].ligand_resname,
                                                        structure = None,
                                                        ligand_resnumber = Ligand[i].res_number)

                #tries to write the ligand pdb
                try:
                    Ligand[i].ligand_filename = Ligand[i].write_PDB(f'{Protein.protein_id}_{Ligand[i].ligand_resname}_lgand{i}.pdb')
                except TypeError as err:
                    raise TypeError(f'{err.args}\ncannot make ligand pdb file for {Ligand[i].ligand_resname}')
                except Exception as err:
                    raise Exception(f'{err} {err.args}\ncannot make ligand pdb file for {Ligand[i].ligand_resname}')
        else:
            Ligand = None

        #use primadorac to get topology and parameter files for any given ligand
        if self.ligand_elaboration_program == 'primadorac':
            if Ligand != None:
        
                #runs primadorac and returns the Ligand list updated with the prm, tpg and new pdb files
                Ligand = funcs4primadorac.run_primadorac(ligand_list = Ligand,
                                                        primadorac_path = self.ligand_elaboration_program_path,
                                                        ph = self.ph)
                
                #As primadorac renames th ligand resnumber as 1
                #and it could create problems they are renamed as -(i + 1)
                #with maximum number -9
                #Not very robust, for sure not the best solution
                for i, ligand in enumerate(pipeline_functions.get_iterable(Ligand)):
                    with open(ligand.ligand_filename, 'r') as f:
                        lines = f.readlines()
                    
                    with open(ligand.ligand_filename, 'w') as f:
                        for line in lines:
                            line = list(line)
                            line[24] = '-'
                            line[25] = f'{(i + 1) % 10}'
                            line = ''.join(line)
                            
                            f.write(line)

        else:
            raise NotImplementedError(self.ligand_elaboration_program)


        if self.MD_program == None or self.MD_program_path == None:
            raise Exception('Need a MD program and a path')

        elif self.MD_program == 'gromacs':

            #makes the necessary resname substitutions for the ForceField
            Protein = funcs4gromacs.residue_substitution(Protein, 'standard')
            Protein = pipeline_functions.get_seqres_PDB(Protein)

        elif self.MD_program == 'orac':

            #makes the necessary resname substitutions for the ForceField
            Protein = funcs4orac.residue_substitution(Protein, 'custom_zinc')
            Protein = pipeline_functions.get_seqres_PDB(Protein)

            #Creating a joined pdb of protein + ligand
            Protein = funcs4orac.join_ligand_and_protein_pdb(Protein, Ligand)

            #first structure optimization, with standard tpg and prm (inside lib module)
            first_opt = funcs4orac.OracFirstOptimization(output_filename = f'first_opt_{Protein.protein_id}.in',
                                                        Protein = Protein,
                                                        Ligand = Ligand,
                                                        protein_tpg_file = self.protein_tpg_file,
                                                        protein_prm_file = self.protein_prm_file,
                                                        MD_program_path = self.MD_program_path)


            Protein.filename = first_opt.execute()

            if Ligand == None:
                raise TypeError('I could not find organic ligands in the structure\n\
                                Maybe the ones you where looking for are in the important_lists.trash list\n\
                                In any case I did some optimizations on your protein, hoping it will come in handy')
        
            solv_box = funcs4orac.OracSolvBoxInput(output_filename = f"{Protein.protein_id}_solv_box.in",
                                                Protein = Protein,
                                                Ligand = Ligand,
                                                protein_tpg_file = self.protein_tpg_file,
                                                protein_prm_file = self.protein_prm_file,
                                                MD_program_path = self.MD_program_path,
                                                solvent_pdb = self.solvent_pdb)
            
            Protein.filename = solv_box.execute()
        
        else:
            raise NotImplementedError(self.MD_program)
        
        


class Ligand_Pipeline(Pipeline):
    """Protein is given as a mmcif
    The ligand is given as input"""
    pass
