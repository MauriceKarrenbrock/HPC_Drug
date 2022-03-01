######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
Contains the classes and functions to create an Hamiltonian Replica Exchange HREM on an HPC cluster with Gromacs
"""

import math
import shutil
from pathlib import Path

import mdtraj

import FSDAMGromacs.get_pbc_atom as _pbc_atom

import HREMGromacs.solute_tempering as _solute_tempering
import HREMGromacs.mdp_files as _hrem_mdp
import HREMGromacs.multidir_setup as _multidir_setup

import PythonPDBStructures.geometry as _pdb_geo

from HPC_Drug.MD import workload_managers
from HPC_Drug.files_IO import write_on_files
from HPC_Drug import important_lists
import HPC_Drug.MD.gromacs.add_dummy_atom as _dummy_atom


class ComplexHREMInput(object):
    
    def __init__(self,
                pdb_file,
                top_file,
                ligand_resname,
                ligand_resnumber,
                kind_of_processor = 'skylake',
                number_of_cores_per_node = 64,
                use_gpu = 'auto',
                gpus_per_node = 1,
                number_of_replicas = 8,
                batteries = None,
                n_steps=None,
                timestep=None,
                constraints='h-bonds',
                temperature=298.15,
                top_scaling_basis=0.2):

        self.pdb_file = pdb_file

        self.top_file = top_file

        self.ligand_resname = ligand_resname

        self.ligand_resnumber = ligand_resnumber

        #the input directory that will be created
        self.HREM_dir = "complex_HREM"

        self.mdp_file = "complex_HREM.mdp"

        self.output_tpr_file = "HREM.tpr"

        self.kind_of_processor = kind_of_processor
        self.number_of_cores_per_node = number_of_cores_per_node

        #gromacs has various options to use gpu
        #auto (default) that will use all the available ones automatically
        #cpu uses no GPU even if available
        #gpu forces the use of GPU (but in case you want to use a gpu auto would be safer and more robust)
        self.use_gpu = use_gpu.lower().strip()
        if self.use_gpu not in ('auto', 'cpu', 'gpu'):
            raise ValueError(f"{self.use_gpu} is not a valid gpu option, valid options are auto cpu gpu")

        self.gpus_per_node = gpus_per_node

        if batteries is None:
            
            self.BATTERIES = self._get_BATTERIES()

        else:

            self.BATTERIES = batteries

        #the replicas for BATTERY
        self.replicas = number_of_replicas

        if timestep is None:

            self.timestep = 0.002

        else:

            self.timestep = timestep

        if n_steps is None:

            # 32 ns in total
            self.n_steps = math.ceil((32000. / self.timestep) / self.BATTERIES)

        else:

            self.n_steps = int(n_steps)

        self.constraints = constraints
        if self.constraints is None:

            self.constraints = 'h-bonds'

        self.temperature = temperature
        if self.temperature is None:

            self.temperature = 298.15

        self.top_scaling_basis = top_scaling_basis


    def _write_template_on_file(self):
        """Creates the mdp
        
        overwrites superclass method
        """

        protein_pbc_atom = _pbc_atom.get_protein_pbc_atom(self.pdb_file)

        mdp_obj = _hrem_mdp.ProteinLigandHREM(mdp_file=self.mdp_file,
            timestep_ps=self.timestep,
            number_of_steps=self.n_steps,
            temperature=self.temperature,
            COM_pull_goups=[
                'Protein',
                f'{self.ligand_resname}',
                'DUM'
            ],
            harmonic_kappa=[
                [
                    'Protein', f'{self.ligand_resname}', 120
                ],

                [
                    'Protein', 'DUM', 120
                ],

                [
                    f'{self.ligand_resname}', 'DUM', 0
                ]
            ],
            pbc_atoms=(
                protein_pbc_atom,
                0,
                0
            ),
            constraints=self.constraints)

        mdp_obj.execute()


    def _get_ns_per_day(self):
        """
        private
        
        Get's the number of ns per day that a kind_of_processor
        processor can process on the sistem
        """

        number_of_atoms = mdtraj.load(self.pdb_file).n_atoms

        if self.use_gpu == 'cpu':
            ns_per_day = ( 15000. / number_of_atoms ) * important_lists.processor_kind_ns_per_day_15000_atoms_for_cpu_only_runs[self.kind_of_processor]
        
        elif self.use_gpu in ('auto', 'gpu'):
            ns_per_day = ( 15000. / number_of_atoms ) * important_lists.ns_per_day_15000_on_gpu_accellerated_architectures

        else:
            raise ValueError(f"{self.use_gpu} is not a valid gpu option, valid options are auto cpu gpu")

        return ns_per_day

    def _get_BATTERIES(self):
        """Get's the number of batteries for REM"""

        ns_per_day = self._get_ns_per_day()

        BATTERIES = math.ceil( 32. / ns_per_day )

        return BATTERIES


    def _write_workloadmanager_inputs(self):
        """private method called by self.execute()"""

        
        wl_manager = MakeWorkloadManagerInput(tpr_file = self.output_tpr_file,
                                            HREM_dir = self.HREM_dir,
                                            BATTERIES = self.BATTERIES,
                                            replicas = self.replicas,
                                            gpu = self.use_gpu,
                                            gpu_per_node = self.gpus_per_node,
                                            cpus_per_node = self.number_of_cores_per_node)

        wl_manager.execute()

    def _make_TPR_files_script(self, scaled_topologies):
        """
        private
        
        writes a bash script to create the .tpr files in loco on the HPC cluster
        """

        filename = "MAKE_TPR_FILES.sh"

        string = "#!/bin/bash\n \n##THIS SCRIPT CREATES THE TPR FILES RUN IT BEFORE THE WORKLOADMANAGER ONE\n#IT IS FOR A PLUMED PATCHED GROMACS INSTALLATION\n\n\n"

        for i in range(self.BATTERIES):
            for j in range(self.replicas):
                string = string + f"gmx grompp -maxwarn 100 -o BATTERY{i}/scaled{j}/{Path(self.output_tpr_file).name} -f {Path(self.mdp_file).name} -p {Path(scaled_topologies[j]).name} -c {Path(self.pdb_file).name} \n"


        write_on_files.write_file(lines = [string], file_name = self.HREM_dir + "/" + filename)

    def _create_scaled_topologies(self):
        """PRIVATE

        Returns
        ---------
        list(pathlib.Path)
            the scaled topologies
        """

        mdtraj_trajectory = mdtraj.load(str(self.pdb_file))

        resSeq_to_scale = _pdb_geo.get_nearest_neighbors_residues_with_mdtraj(
            mdtraj_trajectory,
            ligand_atoms=f'resSeq {self.ligand_resnumber}')
        resSeq_to_scale = resSeq_to_scale[0] + resSeq_to_scale[1]

        # From resid (0 indexed) to resSeq
        top = mdtraj_trajectory.topology
        resSeq_to_scale = [top.atom(top.select(f'resid {resid}')[0]).residue.resSeq for resid in resSeq_to_scale]

        return _solute_tempering.prepare_topologies_for_hrem(top_file=self.top_file,
                                protein_resSeq_to_scale=resSeq_to_scale,
                                ligand_resname=self.ligand_resname,
                                mdp_file=self.mdp_file,
                                gro_file=self.pdb_file,
                                number_of_replicas=self.replicas,
                                basis=self.top_scaling_basis,
                                plumed_path='plumed',
                                preprocess_top=False)


    def execute(self):
        """This method creates the input to make a REM simulation on a HPC cluster
        
        Returns
        ------------
        self.pdb_file, self.top_file, self.HREM_dir : str, str, str
        """

        mdtraj.load(self.pdb_file).save(self.pdb_file + '.gro')

        self.pdb_file, self.top_file = _dummy_atom.add_heavy_dummy_atom(
            gro_file=self.pdb_file + '.gro',
            top_file=self.top_file,
            output_gro='complex_with_dummy.gro',
            output_top='complex_with_dummy.top')

        Path(self.HREM_dir).mkdir(parents=True, exist_ok=True)

        _multidir_setup.make_multiple_hrem_batteries(number_of_batteries=self.BATTERIES,
            replicas_per_battery=self.replicas,
            plumed_file='empty_plumed.dat',
            directory=self.HREM_dir)

        #write the mdp file
        self._write_template_on_file()
        #and copy it in the HREM dir
        shutil.copy(self.mdp_file, self.HREM_dir)

        scaled_topologies = self._create_scaled_topologies()

        # Parmed continuously renames Protein as system1
        # and when
        # defining a restraint gromacs will complain
        with open(self.top_file, 'r') as f:
            protein_top_str = f.read()

        with open(self.top_file, 'w') as f:
            f.write(protein_top_str.replace('system1', 'Protein'))

        #copy the needed protein files in the new directories
        for j in [self.pdb_file, self.top_file] + scaled_topologies:

            shutil.copy(str(j), self.HREM_dir)

        #write workload manager input for different workload managers (slurm pbs ...)
        #already in the self.HREM_dir
        self._write_workloadmanager_inputs()


        #make and copy the script that will make the tpr files in loco
        self._make_TPR_files_script(scaled_topologies = [i.name for i in scaled_topologies])


        #create a file with important info for post processing of the HREM output
        # key = value
        important_info = [
            f"ligand_resname = {self.ligand_resname}\n", #only takes the first ligand, because actually there should be ony one
            f"ligand_resnum = {self.ligand_resnumber}\n",
            f"top_file = {Path(self.top_file).name}\n"
            f"gro_file = {Path(self.self.pdb_file).name}\n"
        ]

        write_on_files.write_file(lines = important_info, file_name = self.HREM_dir + "/" + "important_info.dat")

        return self.pdb_file, self.top_file, self.HREM_dir


class LigandHREMInput(object):
    """
    Makes the HREM to a system where there is only the ligand
    for any ligand in Protein._ligands

    Solvent can or cannot be present (if you want to use the output for FS-DAM there should be no solvent)
    """

    def __init__(self,
                ligand_gro_file,
                ligand_top_file,
                water_gro_file,
                water_ligand_top_file,
                ligand_resname,
                kind_of_processor = 'skylake',
                number_of_cores_per_node = 64,
                use_gpu = 'auto',
                gpus_per_node = 1,
                number_of_replicas = 8,
                batteries = None,
                n_steps=None,
                timestep=None,
                constraints='h-bonds',
                temperature=298.15,
                top_scaling_basis=0.1):

        self.ligand_gro_file = ligand_gro_file

        self.ligand_top_file = ligand_top_file

        self.water_gro_file = water_gro_file

        self.water_ligand_top_file = water_ligand_top_file

        self.ligand_resname = ligand_resname

        #the input directory that will be created
        self.HREM_dir = "ligand_HREM"

        self.mdp_file = "ligand_HREM.mdp"

        self.output_tpr_file = "HREM.tpr"

        self.kind_of_processor = kind_of_processor
        self.number_of_cores_per_node = number_of_cores_per_node

        #gromacs has various options to use gpu
        #auto (default) that will use all the available ones automatically
        #cpu uses no GPU even if available
        #gpu forces the use of GPU (but in case you want to use a gpu auto would be safer and more robust)
        self.use_gpu = use_gpu.lower().strip()
        if self.use_gpu not in ('auto', 'cpu', 'gpu'):
            raise ValueError(f"{self.use_gpu} is not a valid gpu option, valid options are auto cpu gpu")

        self.gpus_per_node = gpus_per_node

        if batteries is None:
        
            self.BATTERIES = 1

        else:

            self.BATTERIES = batteries

        #the replicas for BATTERY
        self.replicas = number_of_replicas

        if timestep is None:

            self.timestep = 0.002

        else:

            self.timestep = timestep

        if n_steps is None:

            # 32 ns in total
            self.n_steps = math.ceil((32000. / self.timestep) / self.BATTERIES)

        else:

            self.n_steps = int(n_steps)

        if constraints is None:

            constraints = 'h-bonds'

        self.constraints = constraints

        self.temperature = temperature
        if self.temperature is None:

            self.temperature = 298.15

        self.top_scaling_basis = top_scaling_basis

    
    def _write_template_on_file(self):
        """Creates the mdp
        
        overwrites superclass method
        """

        mdp_obj = _hrem_mdp.OnlyLigandHREM(mdp_file=self.mdp_file,
            timestep_ps=self.timestep,
            number_of_steps=self.n_steps,
            temperature=self.temperature,
            COM_pull_goups=None,
            harmonic_kappa=None,
            pbc_atoms=None,
            constraints=self.constraints)

        mdp_obj.execute()


    def _make_TPR_files_script(self, scaled_topologies):
        """
        private
        
        writes a bash script to create the .tpr files in loco on the HPC cluster
        """

        filename = "MAKE_TPR_FILES.sh"

        string = "#!/bin/bash\n \n##THIS SCRIPT CREATES THE TPR FILES RUN IT BEFORE THE WORKLOADMANAGER ONE\n#IT IS FOR A PLUMED PATCHED GROMACS INSTALLATION\n\n\n"

        for i in range(self.BATTERIES):
            for j in range(self.replicas):
                string = string + f"gmx grompp -maxwarn 100 -o BATTERY{i}/scaled{j}/{Path(self.output_tpr_file).name} -f {Path(self.mdp_file).name} -p {Path(scaled_topologies[j]).name} -c {Path(self.ligand_gro_file).name} \n"


        write_on_files.write_file(lines = [string], file_name = self.HREM_dir + "/" + filename)


    def _create_scaled_topologies(self):
        """PRIVATE

        Returns
        ---------
        list(pathlib.Path)
            the scaled topologies
        """
        mdtraj_trajectory = mdtraj.load(str(self.ligand_gro_file))
        top = mdtraj_trajectory.topology
        resSeq_to_scale = [top.atom(top.select('all')[0]).residue.resSeq]

        return _solute_tempering.prepare_topologies_for_hrem(top_file=self.ligand_top_file,
                                protein_resSeq_to_scale=resSeq_to_scale,
                                ligand_resname=self.ligand_resname,
                                mdp_file=self.mdp_file,
                                gro_file=self.ligand_gro_file,
                                number_of_replicas=self.replicas,
                                basis=self.top_scaling_basis,
                                plumed_path='plumed',
                                preprocess_top=False)


    def _write_workloadmanager_inputs(self):
        """private method called by self.execute()"""

        
        wl_manager = MakeWorkloadManagerInput(tpr_file = self.output_tpr_file,
                                            HREM_dir = self.HREM_dir,
                                            BATTERIES = self.BATTERIES,
                                            replicas = self.replicas,
                                            gpu = self.use_gpu,
                                            gpu_per_node = self.gpus_per_node,
                                            cpus_per_node = self.number_of_cores_per_node)

        wl_manager.execute()


    def execute(self):
        """This method creates the input to make a REM simulation on a
        HPC cluster
        
        Returns
        ---------
        self.HREM_dir : str
        """

        #the input directory that will be created
        Path(self.HREM_dir).mkdir(parents=True, exist_ok=True)

        _multidir_setup.make_multiple_hrem_batteries(number_of_batteries=self.BATTERIES,
            replicas_per_battery=self.replicas,
            plumed_file='empty_plumed.dat',
            directory=self.HREM_dir)

        #write the mdp file
        self._write_template_on_file()
        #and copy it in the HREM dir
        shutil.copy(self.mdp_file, self.HREM_dir)

        scaled_topologies = self._create_scaled_topologies()

        #copy the needed files in the new directories
        for j in [self.ligand_gro_file,
            self.ligand_top_file,
            self.water_ligand_top_file,
            self.water_gro_file] + scaled_topologies:

            shutil.copy(str(j), self.HREM_dir)
        
        #write workload manager input for different workload managers (slurm pbs ...)
        #already in the self.HREM_dir
        self._write_workloadmanager_inputs()
        
        #make and copy the script that will make the tpr files in loco
        self._make_TPR_files_script(scaled_topologies=scaled_topologies)

        #create a file with important info for post processing of the HREM output
        # key = value
        important_info = [
            f"ligand_resname = {self.ligand_resname}\n",
            f"top_file = {Path(self.ligand_top_file).name}\n",
            f"only_solvent_gro = {Path(self.water_gro_file).name}\n",
            f"solvated_top_file = {Path(self.water_ligand_top_file).name}\n"
            f"gro_file = {Path(self.ligand_gro_file).name}\n"
        ]

        write_on_files.write_file(lines=important_info, file_name=self.HREM_dir + "/" + "important_info.dat")

        return self.HREM_dir





class MakeWorkloadManagerInput(object):

    def __init__(self,
                tpr_file,
                HREM_dir,
                BATTERIES = 1,
                replicas = 8,
                gpu = "auto",
                gpu_per_node = 1,
                cpus_per_node = 64):


        self.tpr_file = tpr_file
        self.HREM_dir = HREM_dir
        self.BATTERIES = BATTERIES
        self.replicas = replicas
        self.gpu = gpu
        self.gpu_per_node = gpu_per_node
        self.cpus_per_node = cpus_per_node


    def write_gromacs_mpirun_string(self, cpus_per_task):
        """private"""

        def gpu_options(use_gpu):
            if use_gpu == 'auto':
                return ""

            else:
                string = f" -nb {use_gpu} -pme {use_gpu} -bonded {use_gpu} -update {use_gpu} -pmefft {use_gpu}"
                
                if use_gpu == 'gpu':
                    string = string + " -npme 1"

                return string

        def multidir_string(battery):
            
            string = " -multidir "

            for i in range(self.replicas):
                string = string + f" BATTERY{battery}/scaled{i} "

            return string


        string = ""            

        for i in range(self.BATTERIES):
            
            string = string + f"mpirun -np {cpus_per_task * self.replicas} gmx_mpi mdrun {gpu_options(use_gpu = self.gpu)} -v -plumed empty_plumed.dat -replex 100 -hrex -dlb no {multidir_string(i)} -s {self.tpr_file} -deffnm HREM & \n"

        string = string + "wait\n"

        return string

    def execute(self):

        if self.gpu == "cpu":

            cpus_per_task = 8
            tasks = self.BATTERIES * self.replicas
            nodes = math.ceil((tasks) / (math.floor((self.cpus_per_node) / (cpus_per_task))))
            wall_time = "24:00:00"
            tasks_per_node = math.floor(self.cpus_per_node / cpus_per_task)
            GPUs = None
            output = "HREM_stdout.out"
            error = "HREM_stderr.err"

            mpirun_string = self.write_gromacs_mpirun_string(cpus_per_task = cpus_per_task)

        else:

            tasks = self.BATTERIES * self.replicas
            nodes = math.ceil((tasks) / self.gpu_per_node)
            wall_time = "24:00:00"
            tasks_per_node = self.gpu_per_node
            cpus_per_task = math.ceil(self.cpus_per_node / tasks_per_node)
            GPUs = tasks
            output = "HREM_stdout.out"
            error = "HREM_stderr.err"

            mpirun_string = self.write_gromacs_mpirun_string(cpus_per_task = cpus_per_task)

        slurm = workload_managers.SlurmHeader(nodes = nodes,
                                            tasks = tasks,
                                            tasks_per_node = tasks_per_node,
                                            cpus_per_task = cpus_per_task,
                                            GPUs = GPUs,
                                            wall_time = wall_time,
                                            output = output,
                                            error = error)

        slurm_header = slurm.execute()

        slurm_header.append(mpirun_string)

        slurm_header = ["\n".join(slurm_header)]

        write_on_files.write_file(lines = slurm_header, file_name = self.HREM_dir + "/" + f"HREM_input.slr")



        pbs = workload_managers.PBSHeader(nodes = nodes,
                                            tasks = tasks,
                                            tasks_per_node = tasks_per_node,
                                            cpus_per_task = cpus_per_task,
                                            GPUs = GPUs,
                                            wall_time = wall_time,
                                            output = output,
                                            error = error)

        pbs_header = pbs.execute()

        pbs_header.append(mpirun_string)

        pbs_header = ["\n".join(pbs_header)]

        write_on_files.write_file(lines = pbs_header, file_name = self.HREM_dir + "/" + f"HREM_input.pbs")
