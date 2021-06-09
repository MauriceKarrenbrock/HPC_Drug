######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

import os
import os.path
import subprocess
import copy
import shutil
import math
import random

import mdtraj

import PythonPDBStructures.trajectories.extract_frames as _extract_frames

import FSDAMGromacs.get_pbc_atom as _pbc_atom
import FSDAMGromacs.pipelines.preprocessing as _preprocessing

from HPC_Drug.auxiliary_functions import path
from HPC_Drug.files_IO import read_file, write_on_files
from HPC_Drug.MD import workload_managers
from HPC_Drug.MD.gromacs import gro2pdb
from HPC_Drug.PDB import biopython, prody
from HPC_Drug import orient


class FSDAMInputPreprocessing(object):

    def __init__(self,
                gromacs_path = "gmx",
                vdw_timestep_ps=None,
                q_timestep_ps=None,
                vdw_number_of_steps=None,
                q_number_of_steps=None,
                creation=True,
                pbc_atoms=None,
                number_of_frames_to_use=200,
                constrains=None):

        #the root directory of the HREM
        self.HREM_dir = os.getcwd()

        self.gromacs_path = gromacs_path

        if vdw_timestep_ps is not None:
            vdw_timestep_ps = float(vdw_timestep_ps)
        self.vdw_timestep_ps = vdw_timestep_ps

        if q_timestep_ps is not None:
            q_timestep_ps = float(q_timestep_ps)
        self.q_timestep_ps = q_timestep_ps

        if vdw_number_of_steps is not None:
            vdw_number_of_steps = int(vdw_number_of_steps)
        self.vdw_number_of_steps = vdw_number_of_steps

        if q_number_of_steps is not None:
            q_number_of_steps = int(q_number_of_steps)
        self.q_number_of_steps = q_number_of_steps

        self.creation = creation

        self.pbc_atoms = pbc_atoms

        self.number_of_frames_to_use = int(number_of_frames_to_use)

        self.constrains = constrains


    def _create_restart_configs(self, fsdam_dir):
        """
        Private

        Makes the starting files for fs-dam and writes them in the RESTART dir
        """

        output_files = []
        i = 0
        while os.path.exists(f"BATTERY{i}"):

            number_of_files = _extract_frames.extract_frames(
                delta_steps=100,
                trajectory=f'BATTERY{i}/scaled0/HREM.trr',
                topology=f'BATTERY{i}/scaled0/HREM.tpr',
                output_name=f'{fsdam_dir}/BATTERY{i}_',
                output_format='gro',
                starting=100,
                end=None
            )

            for k in range(number_of_files):

                output_files.append(f'{fsdam_dir}/BATTERY{i}_{k}.gro')
            
            i = i + 1

        return output_files

    def _make_input_files(self, ligand_resname, top_file, starting_structures):

        for i, item in enumerate(starting_structures):

            starting_structures[i] = item.lsplit('/', 1)[-1]


        if self.creation:

            COM_pull_groups=None

            harmonic_kappa = None

        else:

            COM_pull_groups=['Protein', ligand_resname, 'DUM']

            harmonic_kappa=[
                [
                    'Protein', ligand_resname, 120
                ],

                [
                    'Protein', 'DUM', 120
                ],

                [
                    ligand_resname, 'DUM', 0
                ]
            ]


        output_dictionary = _preprocessing.PreprocessGromacsFSDAM(
            topology_files=[top_file],
            md_program_path=self.gromacs_path,
            alchemical_residue=ligand_resname,
            structure_files=starting_structures,
            COM_pull_groups=COM_pull_groups,
            harmonic_kappa=harmonic_kappa,
            temperature=298.15,
            pbc_atoms=self.pbc_atoms,
            creation=self.creation,
            vdw_timestep_ps = self.vdw_timestep_ps,
            q_timestep_ps = self.q_timestep_ps,
            vdw_number_of_steps = self.vdw_number_of_steps,
            q_number_of_steps = self.q_number_of_steps,
            constrains=self.constrains
            )

        return output_dictionary

    def _make_TPR_files_script(self, q_lines, vdw_lines, fsdamdir):

        """
        private
        
        writes a bash script to create the .tpr files in loco on the HPC cluster
        """

        #Q TPR file
        filename = f"{fsdamdir}/MAKE_Q_TPR_FILES.sh"

        string = "#!/bin/bash\n \n##THIS SCRIPT CREATES THE Q TPR FILES RUN IT BEFORE THE WORKLOADMANAGER ONE\n"
        string += "if you are annihilating make first Q then VDW, if you are creating viceversa\n\n\n"

        string += '\n'.join(q_lines)
        write_on_files.write_file(lines = [string], file_name = filename)

        #VDW TPR file
        filename = f"{fsdamdir}/MAKE_VDW_TPR_FILES.sh"

        string = "#!/bin/bash\n \n##THIS SCRIPT CREATES THE VDW TPR FILES RUN IT BEFORE THE WORKLOADMANAGER ONE\n"
        string += "if you are annihilating make first Q then VDW, if you are creating viceversa\n\n\n"

        string += '\n'.join(vdw_lines)
        write_on_files.write_file(lines = [string], file_name = filename)

    def _select_frames_to_use(self, files, not_used_dir):

        #randomly shuffle the structure files
        random.shuffle(files)

        for _ in range(len(files) - self.number_of_frames_to_use):

            shutil.move(str(files.pop(-1)), str(not_used_dir), copy_function=shutil.copy)

        return files
 
    @staticmethod
    def _add_water_box(water_box, structures):

        water = mdtraj.load(water_box)
        water.center_coordinates()
        #water.make_molecules_whole()

        #the ligand is put in a box of solvated waters
        #and it's file overwritten with the new one
        for i in structures:

            ligand = mdtraj.load(str(i))

            ligand.center_coordinates()

            #ligand.make_molecules_whole()

            joined = water.stack(ligand) #box info are taken from left operand

            joined.save(str(i), force_overwrite=True)

    @staticmethod
    def _edit_top(ligand_top, only_solvent_top):

        only_solvent_lines = read_file.read_file(file_name = only_solvent_top)

        for i in range(len(only_solvent_lines) - 1, -1, -1):

            if only_solvent_lines[i].strip():

                solvent_line = only_solvent_lines[i].rstrip()

                break

        del only_solvent_lines

        ligand_lines = read_file.read_file(file_name = ligand_top)

        for i in range(len(ligand_lines) - 1, -1, -1):

            if ligand_lines[i].strip():

                ligand_lines[i] = solvent_line + '\n' + ligand_lines[i]

                break

        write_on_files.write_file(lines = ligand_lines, file_name = ligand_top)


    def execute(self):

        if self.HREM_dir != os.getcwd():

            os.chdir(self.HREM_dir)

        fsdam_dir = "RESTART"

        os.makedirs(fsdam_dir, exist_ok=True)
        os.makedirs(f'{fsdam_dir}/not_used_frames', exist_ok=True)

        #making some defaults
        useful_info = {
            "ligand_resname" : 'LIG',
            "top_file" : "topol.top",
            "ligand_itp" : "LIG.itp",
            f"only_solvent_gro" : None,
            f"only_solvent_top" : None
        }

        lines = read_file.read_file(file_name = "important_info.dat")
        for line in lines:

            line = line.strip()

            if line:

                line = line.split("=")
                useful_info[line[0].strip()] = line[1].strip()

        starting_configurations = self._create_restart_configs(fsdam_dir = fsdam_dir)

        starting_configurations = self._select_frames_to_use(starting_configurations, not_used_dir=f'{fsdam_dir}/not_used_frames')

        if self.creation:
            
            self._add_water_box(useful_info['only_solvent_gro'], starting_configurations)

            shutil.copy(useful_info["top_file"], 'ligand_solvent_topology.top')

            useful_info["top_file"] = 'ligand_solvent_topology.top'

            self._edit_top(ligand_top=useful_info["top_file"],
                only_solvent_top=useful_info["only_solvent_top"])

        output_dictionary = self._make_input_files(
            ligand_resname=useful_info["ligand_resname"],
            top_file=useful_info["top_file"],
            starting_structures=starting_configurations,
        )

        #writes the TPR scripts in the FSDAM dir
        self._make_TPR_files_script(
            q_lines=output_dictionary['make_q_tpr'],
            vdw_lines=output_dictionary['make_vdw_tpr'],
            fsdamdir=fsdam_dir
        )

        #lazily copy all the mdp, itp and top files in the new dir
        for i in os.listdir(os.getcwd()):

            if i[-4:] in ('.mdp', '.itp', '.top'):

                shutil.copy(i, fsdam_dir)


        #write a suggestion of run file
        with open(f'{fsdam_dir}/RUN_Q.sh', 'w') as f:

            number_of_runs = len(output_dictionary['run_q'])
            f.write(f'# there are {number_of_runs} runs to do\n')
            f.write('# I suggest to use one GPU per run\n\n\n')

            f.write('\n'.join(output_dictionary['run_q']))


        with open(f'{fsdam_dir}/RUN_VDW.sh', 'w') as f:

            number_of_runs = len(output_dictionary['run_vdw'])
            f.write(f'# there are {number_of_runs} runs to do\n')
            f.write('# I suggest to use one GPU per run or even better a job array\n\n\n')

            f.write('\n'.join(output_dictionary['run_vdw']))