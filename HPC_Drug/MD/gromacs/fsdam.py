######################################################################################
# Copyright (c) 2020-2023 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

import os
import os.path
import shutil
import random
from pathlib import Path
import warnings
from multiprocessing import Pool

import mdtraj

import PythonPDBStructures.trajectories.extract_frames as _extract_frames
import PythonPDBStructures.geometry as _geometry

import FSDAMGromacs.pipelines.preprocessing as _preprocessing
import PythonFSDAM.safety_checks as _safety_checks

from HPC_Drug.files_IO import read_file, write_on_files

def _frame_extract_helper(battery, not_used_dir):
    """helper function for parallelization

    I extracted it from the class to have a little
    problems due to pickling as possible
    """

    number_of_files = _extract_frames.extract_all_frames(
            trajectory=f'{battery}/scaled0/HREM.trr',
            topology=f'{battery}/scaled0/HREM.tpr',
            output_name=f'{not_used_dir}/{battery}_',
            output_format='gro'
        )

    files = []
    for k in range(number_of_files):

        files.append(Path(f'{not_used_dir}/{battery}_{k}.gro').resolve())

    return files


def create_complex_lig_creation_frames(apo_files,
                                    holo_files,
                                    ligand_in_vacuum_files,
                                    apo_files_top, 
                                    holo_files_top,
                                    ligand_in_vacuum_files_top,
                                    lig_selection_string,
                                    output_dir=".",
                                    verbose=False):
    """Create the starting configurations to create a ligand in a apo protein

    This function will create the input files needed to create a ligand in a protein pocket
    it works by superimposing frames from the replica exchange of the ligand
    in vacuum on the binding pocket of the apo protein by using frames of the
    replica exchange of the holo system as a reference
    All molecules must be whole

    Parameters
    --------------
    apo_files : str or iterable(str)
        Files of the apo replica exchange run(s)
        any format supported by mdtraj is accepted (xtc, trr, pdb, gro, ...)
    holo_files : str or iterable(str)
        Files of the holo replica exchange run(s)
        any format supported by mdtraj is accepted (xtc, trr, pdb, gro, ...)
    ligand_in_vacuum_files : str or iterable(str)
        Files of the ligand in vacuum replica exchange run(s)
        any format supported by mdtraj is accepted (xtc, trr, pdb, gro, ...)
    apo_files_top : str
        A file of the apo system to use as topology for mdtraj (like pdb or gro files)
    holo_files_top : str
        A file of the holo system to use as topology for mdtraj (like pdb or gro files)
    ligand_in_vacuum_files_top : str
        A file of the ligand in vacuum system to use as topology for mdtraj (like pdb or gro files)
    lig_selection_string : str
        A mdtraj selection string that allows to select the ligand
    output_dir : str, default=current working directory
        The directory in which the output files should be saved
        (will be created if it does not exist)
    verbose : bool, defaul=False

    Returns
    ----------------
    int
        the number of files created

    Warning
    ----------
    This function is a bit naive and sometimes the ligand does not get put in the binding pocket
    (I do not know why), therefore check the output files carefully
    """

    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    apo_files = mdtraj.load(apo_files, top=apo_files_top)

    holo_files = mdtraj.load(holo_files, top=holo_files_top)

    ligand_in_vacuum_files = mdtraj.load(ligand_in_vacuum_files,
                                        top=ligand_in_vacuum_files_top)

    holo_files.atom_slice(holo_files.top.select("not water"), inplace=True)

    holo_files.center_coordinates()

    apo_files.center_coordinates()

    ligand_in_vacuum_files.center_coordinates()

    if verbose:
        print("Loading of the apo holo and lig files finished")
        print("Proceding with the creation of the output files")

    i=0
    for apo, holo, lig in zip(apo_files, holo_files, ligand_in_vacuum_files):

        holo.unitcell_vectors = apo.unitcell_vectors

        holo.center_coordinates()

        lig.unitcell_vectors = apo.unitcell_vectors

        lig.center_coordinates()

        holo = holo.superpose(apo, atom_indices=holo.top.select("name CA"), ref_atom_indices=apo.top.select("name CA"))

        holo = holo.atom_slice(holo.top.select(lig_selection_string)) # keep only ligand

        lig = lig.superpose(holo, atom_indices=lig.top.select("not (symbol H)"), ref_atom_indices=holo.top.select("not (symbol H)"))

        apo = apo.stack(lig)

        apo.save(str(output_dir / f"starting_config_{i}.pdb"))

        i += 1

    if verbose:
        print(f"Number of files created = {i}")

    return i



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
                constrains=None,
                reference_frame=None,
                atoms_in_pocket=None,
                atoms_in_pocket_tollerance=0,
                extra_frames=False,
                add_water=False,
                kind_of_system='protein-ligand'): # possible options are "protein-ligand", "only-ligand", "solvated-ligand"

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

        if constrains is not None:
            warnings.warn('Constraints is not implemented yet, the default (all-atoms) will be used '
            'instead, you can change the mdp file manually')
        self.constrains = constrains

        self.reference_frame = reference_frame
        self.atoms_in_pocket = atoms_in_pocket

        self.atoms_in_pocket_tollerance = int(atoms_in_pocket_tollerance)

        self.extra_frames = extra_frames

        self.add_water = add_water

        self.kind_of_system = kind_of_system


    def _create_restart_configs(self, fsdam_dir, not_used_dir, ligand_resname):
        """
        Private

        Makes the starting files for fs-dam and writes them in the RESTART dir
        """

        BATTERY_dirs = [str(i) for i in Path('.').glob('BATTERY*')]

        # Extract files from trajectory only the first time
        if not self.extra_frames:

            if not BATTERY_dirs:
                raise RuntimeError("Did not find BATTERY* directories containing the HREM")

            parallel_input_list = [str(not_used_dir)] * len(BATTERY_dirs)
            parallel_input_list = list(zip(BATTERY_dirs, parallel_input_list))

            number_of_cpu = int(os.environ.get('OMP_NUM_THREADS', 1))
            number_of_cpu = min(number_of_cpu, len(BATTERY_dirs))

            with Pool(number_of_cpu) as pool:
                pool.starmap(
                    _frame_extract_helper, parallel_input_list, chunksize=1)

        # If I am creating extra frames and there are no BATERY dirs it will still work
        # Even though some BATTERY might not be used
        # Knowing the BATTERY will give a more uniform sampling
        if BATTERY_dirs:
            starting_files = []

            for i in range(len(BATTERY_dirs)):
                starting_files.append([str(i) for i in Path('.').glob(f'{not_used_dir}/BATTERY{i}_*.gro')])

        else:
            starting_files = [
                [str(i) for i in Path('.').glob(f'{not_used_dir}/*.gro')]
            ]

        # If it is a protein ligand system need to check
        # if the ligand is in the pocket
        if self.kind_of_system == 'protein-ligand':
            out_of_pocket_dir = (Path(fsdam_dir)) / 'ligand_out_of_pocket'
            out_of_pocket_dir.mkdir(parents=True, exist_ok=True)

            if self.reference_frame is None:
                reference_frame = mdtraj.load('BATTERY0/scaled0/HREM.trr',
                    top=str(starting_files[0][0])).slice(0)
            else:
                reference_frame = mdtraj.load(self.reference_frame)

            nearest_neighbors, _ = _geometry.get_nearest_neighbors_residues_with_mdtraj(
            reference_frame,
            ligand_atoms=f'resname {ligand_resname}'
            )

            nearest_neighbors = [str(i) for i in nearest_neighbors]
            nearest_neighbors = 'resid ' + (' '.join(nearest_neighbors))

            if self.atoms_in_pocket is None:
                atoms_in_pocket = _safety_checks.get_atoms_in_pocket(ligand=f'resname {ligand_resname}',
                    pocket=nearest_neighbors,
                    pdb_file=reference_frame)
            else:
                atoms_in_pocket = self.atoms_in_pocket
        
        else:
            out_of_pocket_dir = None


        # I had to do it that messy to save computing time
        # if it is a protein ligand system check as little frames as possible
        # to get the needed number
        # otherwise do not check and take the needed
        # number of frames
        output_files = []
        n_frames_accepted = 0
        while n_frames_accepted < self.number_of_frames_to_use:
            for i in range(len(starting_files)):

                # Randomly shuffle the structure files every time
                # I do it to avoid possible biases in the random shuffling
                starting_files[i] = random.sample(starting_files[i], len(starting_files[i]))

                file_name = starting_files[i].pop(-1)

                if out_of_pocket_dir is not None:
                    lig_in_pocket = _safety_checks.check_ligand_in_pocket(
                        ligand=f'resname {ligand_resname}',
                        pocket=nearest_neighbors,
                        pdb_file=file_name,
                        n_atoms_inside=atoms_in_pocket,
                        top=None,
                        make_molecules_whole=False,
                        tollerance=self.atoms_in_pocket_tollerance)
                
                else:
                    lig_in_pocket = True

                if lig_in_pocket:
                    file_name = shutil.move(str(file_name), str(fsdam_dir))
                    file_name = Path(file_name).resolve()

                    output_files.append(file_name)

                    n_frames_accepted += 1

                else:
                    shutil.move(str(file_name), str(out_of_pocket_dir))

                if n_frames_accepted == self.number_of_frames_to_use:
                    break

        return output_files


    def _make_input_files(self, ligand_resname, top_file, starting_structures):

        if self.kind_of_system == 'protein-ligand':

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
        
        else:

            COM_pull_groups=None

            harmonic_kappa = None


        fsdam_preprocessing_obj = _preprocessing.PreprocessGromacsFSDAM(
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
            constrains=self.constrains)
            
        output_dictionary = fsdam_preprocessing_obj.execute()

        return output_dictionary

    def _make_TPR_files_script(self, q_lines, vdw_lines, fsdamdir):

        """
        private
        
        writes a bash script to create the .tpr files in loco on the HPC cluster
        """

        #Q TPR file
        filename = f"{fsdamdir}/MAKE_Q_TPR_FILES.sh"

        string = "#!/bin/bash\n \n##THIS SCRIPT CREATES THE Q TPR FILES RUN IT BEFORE THE WORKLOADMANAGER ONE\n"
        string += "#if you are annihilating make first Q then VDW, if you are creating viceversa\n\n\n"

        string += '\n'.join(q_lines)
        string = string.replace(f'{fsdamdir}/', '')
        write_on_files.write_file(lines = [string], file_name = filename)

        #VDW TPR file
        filename = f"{fsdamdir}/MAKE_VDW_TPR_FILES.sh"

        string = "#!/bin/bash\n \n##THIS SCRIPT CREATES THE VDW TPR FILES RUN IT BEFORE THE WORKLOADMANAGER ONE\n"
        string += "#if you are annihilating make first Q then VDW, if you are creating viceversa\n\n\n"

        string += '\n'.join(vdw_lines)
        string = string.replace(f'{fsdamdir}/', '')
        write_on_files.write_file(lines = [string], file_name = filename)

 
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

    def execute(self):

        if self.HREM_dir != os.getcwd():

            os.chdir(self.HREM_dir)

        if self.extra_frames:
            i = 0
            while Path(f'Extra_RESTART{i}').exists():
                i += 1
            Path(f'Extra_RESTART{i}').mkdir()
            fsdam_dir = f'Extra_RESTART{i}'

        else:
            fsdam_dir = "RESTART"
            Path(fsdam_dir).mkdir(exist_ok=True)
            Path(f'{fsdam_dir}/not_used_frames').mkdir(exist_ok=True)
        
        not_used_dir = 'RESTART/not_used_frames'

        #making some defaults
        useful_info = {
            "ligand_resname" : 'LIG',
            "top_file" : "topol.top",
            "ligand_itp" : "LIG.itp",
            f"only_solvent_gro" : None,
            f"solvated_top_file" : None
        }

        lines = read_file.read_file(file_name = "important_info.dat")
        for line in lines:

            line = line.strip()

            if line:

                line = line.split("=")
                useful_info[line[0].strip()] = line[1].strip()

        # Extract frames and check if they are out of pocket
        # if self.kind_of_system == protein-ligand
        starting_configurations = self._create_restart_configs(fsdam_dir = fsdam_dir,
            not_used_dir=not_used_dir,
            ligand_resname=useful_info['ligand_resname'])

        if self.add_water:
            
            # The ligand is added at the end of the water box, therefore the topology file
            # should respect this order
            self._add_water_box(useful_info['only_solvent_gro'], starting_configurations)

            useful_info["top_file"] = useful_info["solvated_top_file"]

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
