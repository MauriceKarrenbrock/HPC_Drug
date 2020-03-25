"""
Copyright (c) 2020-2020 Maurice Karrenbrock

This software is open-source and is distributed under the
GNU Affero General Public License v3 (agpl v3) license

A copy of the license must be included with any copy of the program or part of it
"""

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
                                solvent_pdb = input_dict['solvent_pdb'],
                                residue_substitution = input_dict['residue_substitution'],
                                kind_of_processor = input_dict['kind_of_processor'],
                                number_of_cores_per_node = input_dict['number_of_cores_per_node'],
                                use_gpu = input_dict['use_gpu'])

    else:
        # ligand is given

        return Ligand_Pipeline()
        

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
                solvent_pdb = None,
                residue_substitution = 'standard',
                kind_of_processor = 'skylake',
                number_of_cores_per_node = 64,
                use_gpu = 'auto'):
        
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

        self.residue_substitution = residue_substitution.strip()

        self.kind_of_processor = kind_of_processor
        self.number_of_cores_per_node = number_of_cores_per_node

        #gromacs has various options to use gpu
        #auto (default) that will use all the available ones automaticly
        #cpu uses no GPU even if available
        #gpu forces the use of GPU (but in case you want to use a gpu auto would be safer and more robust)
        self.use_gpu = use_gpu.lower().strip()
        if self.use_gpu not in ('auto', 'cpu', 'gpu'):
            raise ValueError(f"{self.use_gpu} is not a valid gpu option, valid options are auto cpu gpu")

    def download(self):
        import os

        """If requested in input will download pdb file
        If the given local file doesn't exist will download pdb
        otherwise returns the given path without modifying it"""

        #checks if I got a wrong input
        if self.local.lower() != 'no' and self.local.lower() != 'yes':
            raise ValueError(f'local cannot be {self.local} it can only be "yes" "no" or omitted')

        #if local = no or if the given path doesn't exist it will download the file in the given format (mmcif or pdb)
        if self.local == 'no' or self.local == None or (not os.path.exists(self.protein_filename)):
            
            #downloads the pdb or mmcif returning it's name
            return file_manipulation.download_protein_structure(self.protein_id, self.protein_filetype, None)

        else:
            #protein filename points to an existing file
            #nothing is done
            return self.protein_filename

        


class NoLigand_Pipeline(Pipeline):
    """Protein is given as a mmcif or pdb
    The ligand is already inside the protein file"""
    
    def execute(self):
        """The execution of the pipeline"""
        import os

        #If requested in input will download pdb file
        #If the given local file doesn't exist will download pdb
        #otherwise returns the given path without modifying it
        self.protein_filename = self.download()
        

        #Declaring protein instance
        Protein = structures.Protein(protein_id = self.protein_id,
                                    filename = self.protein_filename,
                                    model = self.model,
                                    chain = self.chain,
                                    file_type = self.protein_filetype)
        

        if self.protein_filetype == 'cif':
            #Parses the mmcif for sulf bonds and organic ligands
            #takes information about the residues binding metals
            #
            #parses the seqres from the mmcif file if possible
            #
            #ligand_resnames is a list containing the ligands
            # resnames and resnumbers [[resname, resnumber], [..., ...], ...]
            Protein, ligand_resnames = pipeline_functions.parse_mmcif(Protein)

        elif self.protein_filetype == 'pdb':

            Protein.substitutions_dict = file_manipulation.get_metal_binding_residues_with_no_header(protein_id = Protein.protein_id,
                                                                                                pdb_file = Protein.filename,
                                                                                                mmcif_file = None,
                                                                                                substitutions_dict = {},
                                                                                                protein_chain = Protein.chain,
                                                                                                protein_model = Protein.model)

            Protein.substitutions_dict, Protein.sulf_bonds = file_manipulation.get_disulf_bonds_with_no_header(protein_id = Protein.protein_id,
                                                                                                pdb_file = Protein.filename,
                                                                                                mmcif_file = None,
                                                                                                substitutions_dict = Protein.substitutions_dict,
                                                                                                protein_chain = Protein.chain,
                                                                                                protein_model = Protein.model)

            ligand_resnames = file_manipulation.get_organic_ligands_with_no_header(protein_id = Protein.protein_id,
                                                                                pdb_file = Protein.filename,
                                                                                mmcif_file = None,
                                                                                protein_chain = Protein.chain,
                                                                                protein_model = Protein.model)

        
        repairer = file_manipulation.PDBRepair()
        
        #adds missing atoms and residues
        #changes non standard residues to standard ones
        #takes a PDBx/mmCIF and returns a PDBx/mmCIF
        #Or a PDB and returns a PDB
        Protein.filename = repairer.add_missing_atoms(pdb_id = Protein.protein_id,
                                                    file_type = Protein.file_type,
                                                    input_filename = Protein.filename,
                                                    repairing_method = self.repairing_method,
                                                    output_filename = None,
                                                    ph = self.ph,
                                                    add_H = False)
        
        ########
        #WORK IN PROGRESS
        ############
        #Define the cys_dict (very useful for later when the resnumbers will be messed up)
        #cys_dict = {resnum: cysteine_number}
        # z = file_manipulation.SubstitutionParser()
        # Protein = z.get_cysteine_dict(Protein = Protein)
        # z = None

        #selects only a selected model and chain, and keeps only one conformation for any disordered atom
        Protein = file_manipulation.select_model_chain_custom(Protein = Protein)

        Protein = file_manipulation.mmcif2pdb(Protein = Protein)

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

        #removes the remaining trash ions
        Protein.filename = file_manipulation.remove_trash_metal_ions(pdb_file = Protein.filename)
        
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
                    Ligand[i].ligand_filename = Ligand[i].write_PDB(f'{Ligand[i].ligand_resname}_{Protein.protein_id}_lgand{i}.pdb')
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
            Protein = funcs4gromacs.residue_substitution(Protein, self.residue_substitution)



            #a patch to add possible missing TER lines in pdb files
            #adressing issues #2733 and #2736 on biopython github
            #Will be removed when the issues will be solved
            with open(Protein.filename, 'r') as f:
                lines = f.readlines()

            with open(Protein.filename, 'w') as f:
                for i in range(len(lines)-1):
                    if lines[i][0:4] == 'ATOM' and lines[i + 1][0:6] == 'HETATM':
                        #writes a TER line on the last ATOM of each chain (not robust for not natural AA)
                        lines[i] = lines[i].strip() + '\n' + "{0:<6}{1}{2}{3}".format("TER", lines[i][6:11], ' '*6, lines[i][17:26])
                    
                    f.write(f"{lines[i].strip()}\n")

                f.write(f"{lines[len(lines)-1].strip()}\n")
            #END OF PATCH

            

            #Make the protein's gro and top file
            gro_top_maker = funcs4gromacs.GromacsMakeProteinGroTop(output_filename = "choices.txt",
                                                    Protein = Protein,
                                                    Ligand = Ligand,
                                                    protein_tpg_file = self.protein_tpg_file,
                                                    solvent_model = '7',
                                                    MD_program_path = self.MD_program_path)
            
            Protein = gro_top_maker.execute()

            #merges the protein and the ligand in a single gro and top file
            gro_top_merger = funcs4gromacs.GromacsMakeJoinedProteinLigandTopGro(output_filename = f"{Protein.protein_id}_joined.pdb",
                                                                                Protein = Protein,
                                                                                Ligand = Ligand,
                                                                                protein_tpg_file = self.protein_tpg_file,
                                                                                solvent_model = self.solvent_pdb,
                                                                                MD_program_path = self.MD_program_path)

            Protein = gro_top_merger.execute()

            #makes a fast geometry optimization
            first_opt = funcs4gromacs.GromacsFirstOptimization(input_filename = f'{Protein.protein_id}_joined_optimized.gro',
                                                            output_filename = f"{Protein.protein_id}_first_opt.mdp",
                                                            Protein = Protein,
                                                            Ligand = Ligand,
                                                            protein_tpg_file = self.protein_tpg_file,
                                                            solvent_model = self.solvent_pdb,
                                                            MD_program_path = self.MD_program_path)

            Protein.gro_file = first_opt.execute()

            if Ligand == None:
                raise TypeError('I could not find organic ligands in the structure\n\
                                Maybe the ones you where looking for are in the important_lists.trash list\n\
                                In any case I did some optimizations on your protein, hoping it will come in handy')

            elif len(Ligand) > 1:
                raise ValueError(f"Found more than one Ligand {Ligand}")

            #makes and optimizes a solvent box
            solv_box = funcs4gromacs.GromacsSolvBoxInput(f"{Protein.protein_id}_solv_box.gro",
                                        output_filename = f'{Protein.protein_id}_solv_box.mdp',
                                        Protein = Protein,
                                        Ligand = Ligand,
                                        protein_tpg_file = self.protein_tpg_file,
                                        solvent_model = self.solvent_pdb,
                                        MD_program_path = self.MD_program_path,
                                        box_borders = '0.8')

            Protein.gro_file = solv_box.execute()

            #create the REM input
            rem_input = funcs4gromacs.GromacsREMInput(input_filename = f"{Protein.protein_id}_REM",
                                                    output_filename = f"{Protein.protein_id}_REM.mdp",
                                                    Protein = Protein,
                                                    Ligand = Ligand,
                                                    protein_tpg_file = self.protein_tpg_file,
                                                    MD_program_path = self.MD_program_path,
                                                    solvent_model = self.solvent_pdb,
                                                    kind_of_processor = self.kind_of_processor,
                                                    number_of_cores_per_node = self.number_of_cores_per_node,
                                                    use_gpu = self.use_gpu)

            gromacs_rem_input_file = rem_input.execute()

        elif self.MD_program == 'orac':

            #makes the necessary resname substitutions for the ForceField
            Protein = funcs4orac.residue_substitution(Protein, self.residue_substitution)
            
            #update the seqres with the names needed for Orac
            Protein = funcs4orac.custom_orac_seqres_from_PDB(Protein)

            #Creating a joined pdb of protein + ligand
            Protein = funcs4orac.join_ligand_and_protein_pdb(Protein, Ligand)

            ########
            #WORK IN PROGRESS
            ##########
            #Protein = pipeline_functions.update_sulf_bonds(Protein = Protein)

            #first structure optimization, with standard tpg and prm (inside lib module)
            first_opt = funcs4orac.OracFirstOptimization(output_filename = f'{Protein.protein_id}_first_opt.in',
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

            elif len(Ligand) > 1:
                raise ValueError(f"Found more than one Ligand {Ligand}")
        
            solv_box = funcs4orac.OracSolvBoxInput(output_filename = f"{Protein.protein_id}_solvbox.in",
                                                Protein = Protein,
                                                Ligand = Ligand,
                                                protein_tpg_file = self.protein_tpg_file,
                                                protein_prm_file = self.protein_prm_file,
                                                MD_program_path = self.MD_program_path,
                                                solvent_pdb = self.solvent_pdb)
            
            Protein.filename = solv_box.execute()


            #Create the REM input
            rem_input = funcs4orac.OracREMInput(output_filename = f"{Protein.protein_id}_REM.in",
                                                Protein = Protein,
                                                Ligand = Ligand,
                                                protein_tpg_file = self.protein_tpg_file,
                                                protein_prm_file = self.protein_prm_file,
                                                MD_program_path = self.MD_program_path,
                                                solvent_pdb = self.solvent_pdb,
                                                kind_of_processor = self.kind_of_processor,
                                                number_of_cores_per_node = self.number_of_cores_per_node)

            orac_rem_input_file = rem_input.execute()

            
        else:
            raise NotImplementedError(self.MD_program)
        
        


class Ligand_Pipeline(Pipeline):
    """Protein is given as a mmcif
    The ligand is given as input"""
    pass
