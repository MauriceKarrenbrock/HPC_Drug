######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains the various pipelines
"""

import subprocess
import os

from HPC_Drug.PDB import download_pdb
from HPC_Drug.PDB import mmcif2pdb
from HPC_Drug.PDB import structural_information_and_repair
from HPC_Drug.PDB import prody
from HPC_Drug.PDB import remove_trash_metal_ions
from HPC_Drug.PDB import select_model_chain
from HPC_Drug.PDB import remove_disordered_atoms
from HPC_Drug.PDB.organic_ligand import get_ligand_topology
import HPC_Drug.auxiliary_functions.path as auxiliary_functions_path
from HPC_Drug.auxiliary_functions import get_iterable
from HPC_Drug.structures import ligand
from HPC_Drug.structures import protein
from HPC_Drug.structures import get_ligands
from HPC_Drug.structures import update_ligands
from HPC_Drug import orient


def choose_pipeline(*args, **kwargs):
    """Takes the parsed program input
    and returns the right pipeline class"""
    
    input_dict = kwargs

    # Protein is a pdb file
    if (input_dict['ligand_in_protein'] is None or input_dict['ligand_in_protein'] == 'yes'):
        #ligand will be taken from protein file (at your own risk)
        
        return NoLigandPipeline(protein = input_dict['protein'],
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
                                use_gpu = input_dict['use_gpu'],
                                gpu_per_node = input_dict['gpu_per_node'],
                                number_of_hrem_replicas_per_battery_bound =input_dict['number_of_hrem_replicas_per_battery_bound'],
                                number_of_hrem_replicas_per_battery_unbound =input_dict['number_of_hrem_replicas_per_battery_unbound'],
                                bound_batteries = input_dict['bound_batteries'],
                                unbound_batteries = input_dict['unbound_batteries'],
                                n_steps_bound=input_dict['n_steps_bound'],
                                n_steps_unbound=input_dict['n_steps_unbound'],
                                timestep_bound=input_dict['timestep_bound'],
                                timestep_unbound=input_dict['timestep_unbound'])

    else:
        # ligand is given

        return LigandPipeline()
        

class Pipeline(object):
    """General pipeline class"""
    
    def __init__(self,
                protein = None,
                protein_filetype = 'cif',
                local = 'no',
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
                use_gpu = 'auto',
                gpu_per_node = 1,
                number_of_hrem_replicas_per_battery_bound =8,
                number_of_hrem_replicas_per_battery_unbound =8,
                bound_batteries = None,
                unbound_batteries = None,
                n_steps_bound=None,
                n_steps_unbound=None,
                timestep_bound=None,
                timestep_unbound=None):
        
        self.protein_id = protein
        self.protein_filetype = protein_filetype

        if self.protein_filetype == None:
            self.protein_filetype = 'cif'

        self.local = local.strip().lower()
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
        
        #ligand can be a pdb file or a mmcif file
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
        #auto (default) that will use all the available ones automatically
        #cpu uses no GPU even if available
        #gpu forces the use of GPU (but in case you want to use a gpu auto would be safer and more robust)
        self.use_gpu = use_gpu.lower().strip()
        if self.use_gpu not in ('auto', 'cpu', 'gpu'):
            raise ValueError(f"{self.use_gpu} is not a valid gpu option, valid options are auto cpu gpu")


        self.gpu_per_node = gpu_per_node

        self.number_of_hrem_replicas_per_battery_bound = number_of_hrem_replicas_per_battery_bound
        self.number_of_hrem_replicas_per_battery_unbound = number_of_hrem_replicas_per_battery_unbound
        self.bound_batteries = bound_batteries
        self.unbound_batteries = unbound_batteries

        self.n_steps_bound = n_steps_bound
        self.n_steps_unbound = n_steps_unbound
        self.timestep_bound = timestep_bound
        self.timestep_unbound = timestep_unbound

    def get_protein_file(self):

        """If requested in input will download pdb (or mmCIF) file
        If the given local file doesn't exist will download pdb (or mmCIF)
        otherwise returns the given path without modifying it
        
        all the paths are converted to absolute paths!"""

        #checks if I got a wrong input
        if self.local != 'no' and self.local != 'yes':
            raise ValueError(f'local cannot be {self.local} it can only be "yes" "no" or omitted')

        #if local = no or if the given path doesn't exist it will download the file in the given format (mmcif or pdb)
        if self.local == 'no' or self.local == None:
            
            #downloads the pdb or mmcif returning it's name
            path = download_pdb.download(protein_id = self.protein_id, file_type = self.protein_filetype, pdir = None)

            #the absolute path of the downloaded file is given to self.protein_filename
            self.protein_filename = auxiliary_functions_path.absolute_filepath(path = path)

        elif self.local == 'yes':

            #protein filename points to an existing file (if it doesn't auxiliary_functions_path.absolute_filepath will raise a FileNotFoundError)
            #the path is updated as absolute path
            self.protein_filename = auxiliary_functions_path.absolute_filepath(path =self.protein_filename)

        else:
            
            raise ValueError(f"local can only be 'yes' or 'no' not {self.local}")
        


class GetProteinLigandFilesPipeline(Pipeline):
    """
    A pipeline that returns a clean and repaired "protein and ions" PDB and
    a PDB file for any not trash organic lingand
    starting from both a PDB, an mmCIF file or a protein id

    execute is the method to call
    returns a Protein instance
    """

    def execute(self):
        """
        A pipeline that returns a clean and repaired "protein and ions" PDB and
        a PDB file for any not trash organic lingand
        starting from both a PDB, an mmCIF file or a protein id

        returns a Protein instance
        """

        #If requested in input will download pdb file
        #If the given local file doesn't exist raises FileNotFounfError
        #otherwise updates self.protein_filename with the given path
        #all the paths are converted to absolute paths
        self.get_protein_file()
        

        # creating protein instance
        Protein = protein.Protein(protein_id = self.protein_id,
                                    pdb_file = self.protein_filename,
                                    model = self.model,
                                    chain = self.chain,
                                    file_type = self.protein_filetype,
                                    tpg_file = self.protein_tpg_file,
                                    prm_file = self.protein_prm_file)
        

        #Get Protein.substitutions_dict Protein.sulf_bonds
        #repairs the Protein.pdb_file
        #returns a list containing the resnames and resnumbers of organic ligands
        # [[resname, resnum], [...], ...]
        #if there are none will be None item
        Info_rep = structural_information_and_repair.InfoRepair(Protein = Protein, repairing_method = self.repairing_method)

        Protein, ligand_resnames_resnums = Info_rep.get_info_and_repair()

        #remove still present disordered atoms (if any)
        Protein = remove_disordered_atoms.remove_disordered_atoms(Protein = Protein)

        #selects only a selected model and chain (Protein.model Protein.chain)
        Protein = select_model_chain.select_model_chain(Protein = Protein)

        #if the protein was a mmCIF I convert it to PDB
        Protein = mmcif2pdb.mmcif2pdb(Protein = Protein)

        #create the Ligand instances and add them to Protein._ligands
        Protein = get_ligands.get_ligands(Protein = Protein, ligand_resnames_resnums = ligand_resnames_resnums)

        Protein.update_structure(struct_type = "prody")

        prody_select = prody.ProdySelect(structure = Protein.structure)
        
        #gets the protein's structure from the pdb
        #The only HETATM remaining are the metal ions
        Protein.structure = prody_select.protein_and_ions()

        #Write Protein only pdb
        Protein.write(file_name = f"{Protein.protein_id}_protein.pdb", struct_type = 'prody')

        #removes the remaining trash ions
        Protein = remove_trash_metal_ions.remove_trash_metal_ions(Protein = Protein)


        #quick patch, will do it better
        #The structure is put in the reference system of the
        #inertia tensor
        orient_obj = orient.Orient(Protein = Protein)
        _, _, Rot_matrix = orient_obj.calculate_moment_of_intertia_tensor()
        Protein.structure = orient_obj.base_change_structure()
        Protein.write()

        Ligand = Protein.get_ligand_list()
        for i in range(len(Ligand)):
            Ligand[i].update_structure(struct_type = "biopython")

            Ligand[i].structure = orient_obj.base_change_structure(structure = Ligand[i].structure, rot_matrix = Rot_matrix)

            Ligand[i].write()


        return Protein
        




class NoLigandPipeline(Pipeline):
    """Protein is given as a mmcif or pdb
    The ligand is already inside the protein file"""
    
    def execute(self):
        """The execution of the pipeline"""

        #download the protein file if needed
        #repair the file
        #parse for disulf bonds, metal binding residues
        #creates a file for the protein and ions
        #and a file for any non trash organic ligand
        get_protein_pipeline = GetProteinLigandFilesPipeline(
            protein = self.protein_id,
            protein_filetype = self.protein_filetype,
            local = self.local,
            filepath = self.protein_filename,
            Protein_model = self.model,
            Protein_chain = self.chain,
            ph = self.ph,
            repairing_method = self.repairing_method,
            ligand = self.ligand_filename,
            protein_prm_file = self.protein_prm_file,
            protein_tpg_file = self.protein_tpg_file
        )

        Protein = get_protein_pipeline.execute()
        
        #get .itp .tpg .prm ... files for any organic ligand
        Protein = get_ligand_topology.get_topology(
            Protein = Protein,
            program_path = self.ligand_elaboration_program_path,
            tool = self.ligand_elaboration_program,
            ph = self.ph
        )

        
        if self.MD_program == None or self.MD_program_path == None:
            raise Exception('Need a MD program and a path')

        from HPC_Drug.MD import residue_renaming

        #makes the necessary resname substitutions for the ForceField
        residue_renamer = residue_renaming.ResidueRenamer(Protein = Protein,
                                                        MD_program = self.MD_program,
                                                        substitution = self.residue_substitution,
                                                        ph = self.ph)

        Protein = residue_renamer.execute()

        if self.MD_program == 'gromacs':

            from HPC_Drug.MD.gromacs import make_top, first_opt, solv_box, hrem

            
            #a patch to add possible missing TER lines in pdb files
            #adressing issues #2733 and #2736 on biopython github
            #Will be removed when the issues will be solved
            with open(Protein.pdb_file, 'r') as f:
                lines = f.readlines()

            with open(Protein.pdb_file, 'w') as f:
                for i in range(len(lines)-1):
                    if lines[i][0:4] == 'ATOM' and lines[i + 1][0:6] == 'HETATM':
                        #writes a TER line on the last ATOM of each chain (not robust for not natural AA)
                        lines[i] = lines[i].strip() + '\n' + "{0:<6}{1}{2}{3}".format("TER", lines[i][6:11], ' '*6, lines[i][17:26])
                    
                    f.write(f"{lines[i].strip()}\n")

                f.write(f"{lines[len(lines)-1].strip()}\n")
            #END OF PATCH


            

            #Make the protein's gro and top file
            gro_top_maker = make_top.GromacsMakeProteinGroTop(Protein = Protein,
                                                    solvent_model = self.solvent_pdb,
                                                    MD_program_path = self.MD_program_path)
            
            Protein = gro_top_maker.execute()

            #merges the protein and the ligand in a single gro and top file
            gro_top_merger = make_top.GromacsMakeJoinedProteinLigandTopGro(Protein = Protein,
                                                                        MD_program_path = self.MD_program_path)

            Protein = gro_top_merger.execute()

            #makes a gro and top for the ligands
            Protein = update_ligands.update_ligands(Protein = Protein, chain_model_selection = False)

            only_ligand_top_gro_maker = make_top.GromacsMakeOnlyLigandTopGro(Protein = Protein,
                                                        MD_program_path = self.MD_program_path,
                                                        solvent_model = self.solvent_pdb)

            Protein = only_ligand_top_gro_maker.execute()



            #makes a fast geometry optimization of the protein ligand system
            first_opt_obj = first_opt.GromacsFirstOptimization(Protein = Protein,
                                                        MD_program_path = self.MD_program_path)

            Protein = first_opt_obj.execute()

            #make fast optimization for the only ligands files
            first_opt_only_ligand = first_opt.GromacsFirstOptimizationOnlyLigand(Protein = Protein,
                                                                    MD_program_path = self.MD_program_path)

            Protein = first_opt_only_ligand.execute()


            if Protein.get_ligand_list() == []:
                raise RuntimeError('I could not find organic ligands in the structure\n\
                Maybe the ones you where looking for are in the important_lists.trash list\n\
                In any case I did some optimizations on your protein, hoping it will come in handy')

            elif len(Protein.get_ligand_list()) > 1:
                raise RuntimeError(f"Found more than one Ligand {Protein.get_ligand_list()}")

            #makes and optimizes a solvent box
            solv_box_obj = solv_box.GromacsSolvBoxInput(Protein = Protein,
                                                MD_program_path = self.MD_program_path,
                                                box_borders = '0.8')

            Protein = solv_box_obj.execute()

            #create the REM input for Plumed patched gromacs
            hrem_input = hrem.GromacsHREMInput(Protein = Protein,
                                            MD_program_path = self.MD_program_path,
                                            kind_of_processor = self.kind_of_processor,
                                            number_of_cores_per_node = self.number_of_cores_per_node,
                                            use_gpu = self.use_gpu,
                                            gpus_per_node = self.gpu_per_node,
                                            number_of_replicas = self.number_of_hrem_replicas_per_battery_bound,
                                            batteries = self.bound_batteries,
                                            n_steps=self.n_steps_bound,
                                            timestep=self.timestep_bound)

            Protein = hrem_input.execute()


            #optimize a box of water
            opt_water = solv_box.OptimizeOnlyWaterBox(force_fileld = Protein.tpg_file,
                                                    solvent_model = self.solvent_pdb,
                                                    MD_program_path = self.MD_program_path)

            only_solvent_box_gro, only_solvent_box_top = opt_water.execute()

            #create a HREM dir for the ligand
            ligand_hrem_input = hrem.GromacsHREMOnlyLigand(Protein,
                                                only_solvent_box_gro = only_solvent_box_gro,
                                                only_solvent_box_top = only_solvent_box_top,
                                                MD_program_path = self.MD_program_path,
                                                kind_of_processor = self.kind_of_processor,
                                                number_of_cores_per_node = self.number_of_cores_per_node,
                                                use_gpu = self.use_gpu,
                                                gpus_per_node = self.gpu_per_node,
                                                number_of_replicas = self.number_of_hrem_replicas_per_battery_unbound,
                                                batteries = self.unbound_batteries,
                                                n_steps=self.n_steps_unbound,
                                                timestep=self.timestep_unbound)

            Protein = ligand_hrem_input.execute()

        elif self.MD_program == 'orac':

            from HPC_Drug.MD.orac import first_opt, solv_box, hrem, orac_seqres

            from HPC_Drug.PDB import merge_pdb

            
            #update the seqres with the names needed for Orac
            Protein = orac_seqres.custom_orac_seqres_from_PDB(Protein)

            #Creating a joined pdb of protein + ligand
            Protein = merge_pdb.merge_pdb(Protein = Protein)

            #makes a gro and top for the ligands
            Protein = update_ligands.update_ligands(Protein = Protein, chain_model_selection = False)

            #first structure optimization, with standard tpg and prm (inside lib module)
            first_opt_obj = first_opt.OracFirstOptimization(Protein = Protein,
                                                        MD_program_path = self.MD_program_path)


            Protein = first_opt_obj.execute()

            #make fast optimization for the only ligands files
            first_opt_only_ligand = first_opt.OracFirstOptimizationOnlyLigand(Protein = Protein,
                                                                    MD_program_path = self.MD_program_path)

            Protein = first_opt_only_ligand.execute()

            if Protein.get_ligand_list() == []:
                raise RuntimeError('I could not find organic ligands in the structure\n\
                Maybe the ones you where looking for are in the important_lists.trash list\n\
                In any case I did some optimizations on your protein, hoping it will come in handy')

            elif len(Protein.get_ligand_list()) > 1:
                raise RuntimeError(f"Found more than one Ligand {Protein.get_ligand_list()}")
        
            solv_box_obj = solv_box.OracSolvBoxInput(Protein = Protein,
                                                MD_program_path = self.MD_program_path,
                                                solvent_pdb = self.solvent_pdb)
            
            Protein = solv_box_obj.execute()


            #Create the REM input
            hrem_input_obj = hrem.HREMOracInput(Protein = Protein,
                                                MD_program_path = self.MD_program_path,
                                                solvent_pdb = self.solvent_pdb,
                                                kind_of_processor = self.kind_of_processor,
                                                number_of_cores_per_node = self.number_of_cores_per_node,
                                                number_of_replicas = self.number_of_hrem_replicas_per_battery_bound)

            Protein = hrem_input_obj.execute()

            #optimize a box of only water
            solv_box_only_water = solv_box.OptimizeOnlyWaterBox(solvent_box = None,
                                                                tpg_file = Protein.tpg_file,
                                                                prm_file = Protein.prm_file,
                                                                MD_program_path = self.MD_program_path,
                                                                chain = Protein.chain)

            only_water_pdb = solv_box_only_water.execute()

            #make the hrem input for the ligand
            hrem_only_ligand = hrem.HREMOracInputOnlyLigand(Protein = Protein,
                                                            solvent_box = only_water_pdb,
                                                            MD_program_path = self.MD_program_path,
                                                            number_of_cores_per_node = self.number_of_cores_per_node,
                                                            number_of_replicas = self.number_of_hrem_replicas_per_battery_unbound)

            Protein = hrem_only_ligand.execute()

            
        else:
            raise NotImplementedError(self.MD_program)

        return Protein
        
        


class LigandPipeline(Pipeline):
    """Protein is given as a mmcif
    The ligand is given as input"""
    pass
