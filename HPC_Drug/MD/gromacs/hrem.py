######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
Contains the classes and functions to create an Hamiltonian Replica Exchange HREM on an HPC cluster with Gromacs
"""

import os
import math
import shutil
import subprocess

from HPC_Drug.MD.gromacs import gromacs_input
from HPC_Drug.files_IO import read_file
from HPC_Drug.files_IO import write_on_files
from HPC_Drug.auxiliary_functions import path
from HPC_Drug import orient
from HPC_Drug import important_lists
from HPC_Drug import funcs4slurm
from HPC_Drug import funcs4pbs


class GromacsHREMInput(gromacs_input.GromacsInput):
    #WORK IN PROGRESS
    
    def __init__(self,
                Protein,
                MD_program_path = 'gmx',
                kind_of_processor = 'skylake',
                number_of_cores_per_node = 64,
                use_gpu = 'auto'):

        super().__init__(Protein = Protein,
                        MD_program_path = MD_program_path)


        #the input directory that will be created
        self.HREM_dir = f"{self.Protein.protein_id}_HREM"

        self.BATTERIES = self._get_BATTERIES()
        #the replicas for BATTERY
        self.replicas = 8

        self.elaborated_top_file = f"{self.Protein.protein_id}_elaborated_topology.top"

        self.mdp_file = f"{self.Protein.protein_id}_HREM.mdp"

        self.output_tpr_file = f"{self.Protein.protein_id}_HREM.tpr"

        self.kind_of_processor = kind_of_processor
        self.number_of_cores_per_node = number_of_cores_per_node

        #an instance of orient.Orient class
        self.orient = orient.Orient(self.Protein, self.Protein.get_ligand_list())

        #gromacs has various options to use gpu
        #auto (default) that will use all the available ones automatically
        #cpu uses no GPU even if available
        #gpu forces the use of GPU (but in case you want to use a gpu auto would be safer and more robust)
        self.use_gpu = use_gpu.lower().strip()
        if self.use_gpu not in ('auto', 'cpu', 'gpu'):
            raise ValueError(f"{self.use_gpu} is not a valid gpu option, valid options are auto cpu gpu")


        #as to make a hrem on gromacs you have to trick it in thinking it is doing a temperature rem
        #the temperature will be changed during the process
        self.temperature = 298.15

        self.template = [
            "; VARIOUS PREPROCESSING OPTIONS",
            "; Preprocessor information: use cpp syntax.",
            "; e.g.: -I/home/joe/doe -I/home/mary/roe",
            "include                  =",
            "; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)",
            "define                   =",
            "",
            "; RUN CONTROL PARAMETERS",
            "integrator               = md",
            "; Start time and timestep in ps",

            self._write_TIME_TIMESTEP_string(),

            "; For exact run continuation or redoing part of a run",
            "init-step                = 0",
            "; Part index is updated automatically on checkpointing (keeps files separate)",
            "simulation-part          = 1",
            "; mode for center of mass motion removal",
            "comm-mode                = Linear",
            "; number of steps for center of mass motion removal",
            "nstcomm                  = 100",
            "; group(s) for center of mass motion removal",
            "comm-grps                =",
            "",
            "; TEST PARTICLE INSERTION OPTIONS",
            "rtpi                     = 0.05",
            "",
            "; OUTPUT CONTROL OPTIONS",
            "; Output frequency for coords (x), velocities (v) and forces (f)",
            "nstxout                  = 100000",
            "nstvout                  = 100000",
            "nstfout                  = 100000",
            "; Output frequency for energies to log file and energy file",
            "nstlog                   = 100000",
            "nstcalcenergy            = 100",
            "nstenergy                = 100000",
            "; Output frequency and precision for .xtc file",
            "nstxtcout                = 80000",
            "xtc-precision            = 1000",
            "; This selects the subset of atoms for the .xtc file. You can",
            "; select multiple groups. By default all atoms will be written.",
            "xtc-grps                 =",
            "; Selection of energy groups",
            "energygrps               = System",
            "",
            "; NEIGHBORSEARCHING PARAMETERS",
            "; cut-off scheme (group: using charge groups, Verlet: particle based cut-offs)",
            "; nblist update frequency",
            "cutoff-scheme            = Verlet",
            "nstlist                  = 20",
            "verlet-buffer-tolerance  = 0.0001",
            "; ns algorithm (simple or grid)",
            "ns_type                  = grid",
            "; Periodic boundary conditions: xyz, no, xy",
            "pbc                      = xyz",
            "periodic-molecules       = no",
            "; Allowed energy drift due to the Verlet buffer in kJ/mol/ps per atom,",
            "; a value of -1 means: use rlist",
            "; nblist cut-off",
            "rlist                    = 1.0",
            "; long-range cut-off for switched potentials",
            "rlistlong                = -1",
            "",
            "; OPTIONS FOR ELECTROSTATICS AND VDW",
            "; Method for doing electrostatics",
            "coulombtype              = PME",
            "rcoulomb-switch          = 0",
            "rcoulomb                 = 1.0",
            "; Relative dielectric constant for the medium and the reaction field",
            "epsilon-r                = 1",
            "epsilon-rf               = 0",
            "; Method for doing Van der Waals",
            "vdw-type                 = Cut-off",
            "; cut-off lengths",
            "rvdw-switch              = 0",
            "rvdw                     = 1.0",
            "; Apply long range dispersion corrections for Energy and Pressure",
            "DispCorr                 = EnerPres",
            "; Extension of the potential lookup tables beyond the cut-off",
            "table-extension          = 1",
            "; Separate tables between energy group pairs",
            "energygrp-table          =",
            "; Spacing for the PME/PPPM FFT grid",
            "fourierspacing           = 0.12",
            "; FFT grid size, when a value is 0 fourierspacing will be used",
            "fourier-nx               = 0",
            "fourier-ny               = 0",
            "fourier-nz               = 0",
            "; EWALD/PME/PPPM parameters",
            "pme-order                = 4",
            "ewald-rtol               = 1e-06",
            "ewald-geometry           = 3d",
            "epsilon-surface          =",
            "optimize-fft             = no",
            "",
            "; IMPLICIT SOLVENT ALGORITHM",
            "implicit-solvent         = No",
            "",
            "; OPTIONS FOR WEAK COUPLING ALGORITHMS",
            "; Temperature coupling",
            "tcoupl                   = v-rescale",
            "nsttcouple               = -1",
            "nh-chain-length          = 1",
            "; Groups to couple separately",
            "tc-grps                  = System",
            "; Time constant (ps) and reference temperature (K)",
            "tau-t                    = 0.2",

            f"ref-t                    = {self.temperature}",

            "; pressure coupling",
            "pcoupl                   = Berendsen",
            "pcoupltype               = Isotropic",
            "nstpcouple               = -1",
            "; Time constant (ps), compressibility (1/bar) and reference P (bar)",
            "tau-p                    = 0.5",
            "compressibility          = 4.6e-5",
            "ref-p                    = 1",
            "; Scaling of reference coordinates, No, All or COM",
            "refcoord-scaling         = COM",
            "",
            "; GENERATE VELOCITIES FOR STARTUP RUN",
            "gen-vel                  = no",
            "gen-temp                 = 500",
            "gen-seed                 = 173529",
            "",
            "; OPTIONS FOR BONDS",
            "constraints              = all-bonds",
            "; Type of constraint algorithm",
            "constraint-algorithm     = Lincs",
            "; Do not constrain the start configuration",
            "continuation             = no",
            "; Use successive overrelaxation to reduce the number of shake iterations",
            "Shake-SOR                = no",
            "; Relative tolerance of shake",
            "shake-tol                = 0.00001",
            "; Highest order in the expansion of the constraint coupling matrix",
            "lincs-order              = 5",
            "; Number of iterations in the final step of LINCS. 1 is fine for",
            "; normal simulations, but use 2 to conserve energy in NVE runs.",
            "; For energy minimization with constraints it should be 4 to 8.",
            "lincs-iter               = 2",
            "; Lincs will write a warning to the stderr if in one step a bond",
            "; rotates over more degrees than",
            "lincs-warnangle          = 30",
            "; Convert harmonic bonds to morse potentials",
            "morse                    = no"
        ]


    def _get_ns_per_day(self):
        """
        private
        
        Get's the number of ns per day that a kind_of_processor
        processor can process on the sistem
        """

        number_of_atoms = self.orient.get_first_last_atom_structure(Protein = self.Protein, Ligand = self.Protein.get_ligand_list())
        number_of_atoms = number_of_atoms[0]

        if self.use_gpu == 'cpu':
            ns_per_day = ( 15000. / number_of_atoms ) * important_lists.processor_kind_ns_per_day_15000_atoms_for_cpu_only_runs[self.kind_of_processor]
        
        elif self.use_gpu in ('auto', 'gpu'):
            ns_per_day = ( 15000. / number_of_atoms ) * important_lists.ns_per_day_15000_on_gpu_accellerated_architectures

        else:
            raise ValueError(f"{self.use_gpu} is not a valid gpu option, valid options are auto cpu gpu")

        return ns_per_day

    def _write_TIME_TIMESTEP_string(self):
        """writes the tinit timestep and number of steps string"""

        tinit = 0
        timestep = 0.00150 #ps
        
        #ps calculated in 24 hours by one mpi_run (BATTERIES)
        number_of_steps = self._get_ns_per_day() * 1.E+3
        #number of steps (integer)
        number_of_steps = math.ceil( number_of_steps / timestep )

        string = f"tinit = {tinit} \ndt = {timestep} \nnsteps = {number_of_steps}\n"

        return string

    def _get_BATTERIES(self):
        """Get's the number of batteries for REM"""

        ns_per_day = self._get_ns_per_day()

        BATTERIES = math.ceil( 32. / ns_per_day )

        return BATTERIES

    def _edit_top_file(self, top_file = None):
        """PRIVATE"""

        if top_file is None:
            top_file = self.elaborated_top_file

        #get the hot residues
        hot_residues = self.orient.get_hot_residues_for_rem(Protein = self.Protein, Ligand = self.Protein.get_ligand_list(), cutoff = 4.5, residue_dist = 10.0)
        hot_ids = []
        for residue in hot_residues:
            hot_ids.append(str(residue[1]).strip())

        #get the ligand resnames
        ligands_resnames = []
        for lgand in self.Protein.get_ligand_list():
            ligands_resnames.append(lgand.resname)

        #read the topology
        lines = read_file.read_file(file_name = top_file)


        #heat the right atoms

        #auxiliary bool variable
        is_atoms = False
        for i in range(len(lines)):

            #if the line is empty i can go on with the for loop
            if lines[i].strip() != "":
                tmp_line = lines[i].strip().split()
            else:
                continue

            #check if we are in the atoms part
            if tmp_line[0].replace(" ", "") == "[atoms]":
                is_atoms = True

            #end of atoms part
            elif tmp_line[0][0] == "[":
                is_atoms = False

            if is_atoms:

                #if it is not a comment
                if tmp_line[0][0] != ";":

                    #not heat solvent
                    if tmp_line[3].strip() == "SOL":
                        continue

                    #heat the ligand
                    if tmp_line[3].strip() in ligands_resnames:

                        tmp_line[1] = tmp_line[1] + "_"

                        lines[i] = " ".join(tmp_line)

                    #heat the near residues
                    elif tmp_line[2].strip() in hot_ids:

                        tmp_line[1] = tmp_line[1] + "_"

                        lines[i] = " ".join(tmp_line)

        write_on_files.write_file(lines = lines, file_name = top_file)

    def _plumed_partial_tempering(self, scaling_value, input_file, output_file):
        """
        private

        Creates the scaled topology files with plumed
        """

        plumed = path.absolute_programpath(program = "plumed")

        command = [plumed, "partial_tempering", f"{scaling_value}"]

        print("Running Plumed")

        with open(input_file, "r") as input_file_handle:
            with open(output_file, "w") as output_file_handle:
                
                r = subprocess.run(command,
                                shell = False,
                                stdin = input_file_handle,
                                stdout= output_file_handle,
                                stderr= subprocess.PIPE)

        print(r.stderr)

        if r.returncode != 0:
            raise RuntimeError(f"Plumed failure\n{r.stderr}")

    def _write_workloadmanager_inputs(self):
        """private method called by self.execute()"""

        file_list = []

        #making a nested list the firs dimension reppresents the mpirun
        #and the second the various directories that will be called
        #with each mpirun (gmx -multidir)
        multidir = [[] for x in range(self.BATTERIES)]

        for i in range(self.BATTERIES):
            for j in range(self.replicas):

                multidir[i].append(f"BATTERY{i}/scaled{j}")
 
        slurm = funcs4slurm.SlurmInput(MD_input_file = self.output_tpr_file,
                                    MD_input_directories = multidir,
                                    slurm_input_file = f'{self.mdp_file.rsplit(".", 1)[0]}.slr',
                                    MD_program = 'gromacs',
                                    MD_calculation_type = 'rem',
                                    number_of_cores_per_node = self.number_of_cores_per_node,
                                    max_time = "24:00:00",
                                    ntasks = self.BATTERIES * 8,
                                    cpus_per_task = 8,
                                    std_out = f'{self.mdp_file.rsplit(".", 1)[0]}.out',
                                    std_err = f'{self.mdp_file.rsplit(".", 1)[0]}.err',
                                    use_gpu = self.use_gpu,
                                    plumed = False)

        file_list.append(slurm.write())

        pbs = funcs4pbs.SlurmInput(MD_input_file = self.output_tpr_file,
                                    MD_input_directories = multidir,
                                    slurm_input_file = f'{self.mdp_file.rsplit(".", 1)[0]}.pbs',
                                    MD_program = 'gromacs',
                                    MD_calculation_type = 'rem',
                                    number_of_cores_per_node = self.number_of_cores_per_node,
                                    max_time = "24:00:00",
                                    ntasks = self.BATTERIES * 8,
                                    cpus_per_task = 8,
                                    std_out = f'{self.mdp_file.rsplit(".", 1)[0]}.out',
                                    std_err = f'{self.mdp_file.rsplit(".", 1)[0]}.err',
                                    use_gpu = self.use_gpu,
                                    plumed = False)

        file_list.append(pbs.write())

        return file_list

    def _make_TPR_files_script(self, scaled_topologies):
        """
        private
        
        writes a bash script to create the .tpr files in loco on the HPC cluster
        it makes both the one for a plumed patched gromacs installation and for a non patched one
        """

        # #For NOT plumed patched gromacs

        # filename = "Not_Patched_MAKE_TPR_FILES.sh"

        # string = "#!/bin/bash\n \n##THIS SCRIPT CREATES THE TPR FILES RUN IT BEFORE THE WORKLOADMANAGER ONE\nIT IS FOR A NOT PATCHED GROMACS INSTALLATION\n\n\n"

        # for i in range(self.BATTERIES):
        #     for j in range(self.replicas):
        #         string = string + f"gmx grompp -maxwarn 100 -o BATTERY{i}/scaled{j}/{self.output_tpr_file} -f BATTERY{i}/{self.mdp_file}{j}.mdp -p BATTERY{i}/{scaled_topologies[j]} -c BATTERY{i}/{self.Protein.gro_file.rsplit('/', 1)[-1]} \n"


        # write_on_files.write_file(lines = [string], file_name = self.HREM_dir + "/" + filename)
        

        #For PLUMED PATCHED gromacs

        filename = "Plumed_Patched_MAKE_TPR_FILES.sh"

        string = "#!/bin/bash\n \n##THIS SCRIPT CREATES THE TPR FILES RUN IT BEFORE THE WORKLOADMANAGER ONE\nIT IS FOR A PLUMED PATCHED GROMACS INSTALLATION\n\n\n"

        for i in range(self.BATTERIES):
            for j in range(self.replicas):
                string = string + f"gmx grompp -maxwarn 100 -o BATTERY{i}/scaled{j}/{self.output_tpr_file} -f {self.mdp_file} -p {scaled_topologies[j]} -c {self.Protein.gro_file.rsplit('/', 1)[-1]} \n"


        write_on_files.write_file(lines = [string], file_name = self.HREM_dir + "/" + filename)

        

    def _get_hamiltonian_scaling_values(self):
        """
        Scales the hamiltonian with a geometrical progression
        scale(m) =scale^(m/(nprocs−1)) with scale = 0.2 and nprocs = 8
        0 <= m <= nprocs -1

        for more information check the orac manual http://www.chim.unifi.it/orac/orac-manual.pdf
        page 123 """

        nprocs = 8.0
        scale = 0.2

        #instantiating the list and putting the value for scale^0 = 1.0
        hamiltonian_scaling_values = [1.0]

        for m in range(1, 8):

            scale_m = scale ** ( m / (nprocs -1) )

            hamiltonian_scaling_values.append(scale_m)

        hamiltonian_scaling_values = tuple(hamiltonian_scaling_values)

        return hamiltonian_scaling_values



    def execute(self):
        """This method does not run gromacs but
        creates the input to make a REM simulation on a
        HPC cluster"""

        #creates the REM directory that will be copied to the HPC cluster
        os.makedirs(self.HREM_dir, exist_ok=True)

        #create the BATTERY dirs and the scaled dirs
        for i in range(self.BATTERIES):

            for j in range(self.replicas):

                os.makedirs(self.HREM_dir + "/" + f"BATTERY{i}" + "/" + f"scaled{j}", exist_ok=True)

        #write the mdp file
        self._write_template_on_file()
        #and copy it in the HREM dir
        shutil.copy(self.mdp_file, self.HREM_dir)

        #create the elaborated top file for plumed
        self._interact_with_gromacs(string = [f"{self.MD_program_path}", "grompp", "-f", f"{self.mdp_file}", "-c", f"{self.Protein.gro_file}", "-p", f"{self.Protein.top_file}", "-maxwarn", "100", "-pp", f"{self.elaborated_top_file}"])

        #finds the hot residues and edits the elaborated top file accordingly
        self._edit_top_file(top_file = self.elaborated_top_file)


        #use plumed to make the scaled topologies

        hamiltonian_scaling_values = self._get_hamiltonian_scaling_values()
        scaled_topologies = []
        for counter, value in enumerate(hamiltonian_scaling_values):
            #they are already saved in the HREM dir
            self._plumed_partial_tempering(scaling_value = value,
                                        input_file = self.elaborated_top_file,
                                        output_file = self.HREM_dir + f"{self.Protein.protein_id}_scaled_{counter}.top")

            scaled_topologies.append(f"{self.Protein.protein_id}_scaled_{counter}.top")


        #copy the needed protein files in the new directories
        for j in (self.Protein.gro_file, self.Protein.top_file):

            shutil.copy(j, self.HREM_dir)

        #copy the ligand itp files
        for lgand in self.Protein.get_ligand_list():

            shutil.copy(lgand.itp_file, self.HREM_dir)

        
        #make and copy the workload manager input files
        #write workload manager input for different workload managers (slurm pbs ...)
        workload_files = self._write_workloadmanager_inputs()
        for wl_file in workload_files:
            shutil.copy(wl_file, self.HREM_dir)


        #make and copy the script that will make the tpr files in loco
        self._make_TPR_files_script(scaled_topologies = scaled_topologies)


        #create a file with important info for post processing of the HREM output
        # key = value
        important_info = [
            f"ligand_resname = {self.Protein.get_ligand_list()[0].resname}", #only takes the first ligand, because actually there should be ony one
            f"ligand_itp = {self.Protein.get_ligand_list()[0].itp_file}"
            f"protein_topology = {self.Protein.top_file}"
        ]

        write_on_files.write_file(lines = important_info, file_name = self.HREM_dir + "/" + "important_info.dat")

        return self.Protein



class GromacsHREMOnlyLigand(GromacsHREMInput):
    """
    Makes the HREM to a system where there is only the ligand

    the ligand shall be "delivered" as a Protein instance with an empty _ligands

    Solvent can or cannot be present (if you want to use the output for FS-DAM there should be no solvent)
    """

    def __init__(self,
                Protein,
                MD_program_path = 'gmx',
                kind_of_processor = 'skylake',
                number_of_cores_per_node = 64,
                use_gpu = 'auto'):

        super().__init__(Protein = Protein,
                        MD_program_path = MD_program_path,
                        kind_of_processor = kind_of_processor,
                        number_of_cores_per_node = number_of_cores_per_node,
                        use_gpu = use_gpu)


        #the input directory that will be created
        self.HREM_dir = f"{self.Protein.protein_id}_only_ligand_HREM"

        self.BATTERIES = self._get_BATTERIES()
        #the replicas for BATTERY
        self.replicas = 8

        self.elaborated_top_file = f"{self.Protein.protein_id}_only_ligand_elaborated_topology.top"

        self.mdp_file = f"{self.Protein.protein_id}_only_ligand_HREM.mdp"

        self.output_tpr_file = f"{self.Protein.protein_id}_only_ligand_HREM.tpr"

    def _edit_top_file(self, top_file):

        #read the topology
        lines = read_file.read_file(file_name = top_file)


        #heat the right atoms

        #auxiliary bool variable
        is_atoms = False
        for i in range(len(lines)):

            #if the line is empty i can go on with the for loop
            if lines[i].strip() != "":
                tmp_line = lines[i].strip().split()
            else:
                continue

            #check if we are in the atoms part
            if tmp_line[0].replace(" ", "") == "[atoms]":
                is_atoms = True

            #end of atoms part
            elif tmp_line[0][0] == "[":
                is_atoms = False

            if is_atoms:

                #if it is not a comment
                if tmp_line[0][0] != ";":

                    #not heat solvent
                    if tmp_line[3].strip() != "SOL":
                        tmp_line[1] = tmp_line[1] + "_"

                        lines[i] = " ".join(tmp_line)

        write_on_files.write_file(lines = lines, file_name = top_file)

    def _get_BATTERIES(self):

        #for the ligand one battery is more than enough
        return 1

