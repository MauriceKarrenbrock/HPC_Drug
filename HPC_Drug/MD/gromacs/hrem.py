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

from HPC_Drug.MD import workload_managers
from HPC_Drug.MD.gromacs import gromacs_input
from HPC_Drug.files_IO import read_file
from HPC_Drug.files_IO import write_on_files
from HPC_Drug.auxiliary_functions import path
from HPC_Drug import orient
from HPC_Drug import important_lists


class GromacsHREMInput(gromacs_input.GromacsInput):
    
    def __init__(self,
                Protein,
                MD_program_path = 'gmx',
                kind_of_processor = 'skylake',
                number_of_cores_per_node = 64,
                use_gpu = 'auto',
                gpus_per_node = 1,
                number_of_replicas = 8):

        super().__init__(Protein = Protein,
                        MD_program_path = MD_program_path)


        #the input directory that will be created
        self.HREM_dir = f"{self.Protein.protein_id}_HREM"

        self.elaborated_top_file = f"{self.Protein.protein_id}_elaborated_topology.top"

        self.mdp_file = f"{self.Protein.protein_id}_HREM.mdp"

        self.output_tpr_file = f"HREM.tpr"

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

        self.gpus_per_node = gpus_per_node

        self.BATTERIES = self._get_BATTERIES()
        #the replicas for BATTERY
        self.replicas = number_of_replicas

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

            #if the line is empty or a comment i can go on with the for loop
            if lines[i].strip() != "":
                if lines[i].strip()[0] != ";":
                
                    tmp_line = lines[i].strip().split()

                else:
                    continue
            
            else:
                continue

            #check if we are in the atoms part
            if lines[i].strip().replace(" ", "").split(";")[0] == "[atoms]":
                is_atoms = True

                continue

            #end of atoms part
            elif tmp_line[0][0] == "[":
                is_atoms = False

                continue

            if is_atoms:

                if len(tmp_line) >= 4:

                    #do not heat solvent
                    if tmp_line[3].strip() == "SOL":
                        continue

                    #heat the ligand
                    if tmp_line[3].strip() in ligands_resnames:

                        tmp_line[1] = tmp_line[1] + "_"

                        lines[i] = " ".join(tmp_line) + "\n"

                    #heat the near residues
                    elif tmp_line[2].strip() in hot_ids:

                        tmp_line[1] = tmp_line[1] + "_"

                        lines[i] = " ".join(tmp_line) + "\n"


        write_on_files.write_file(lines = lines, file_name = top_file)

    def _plumed_partial_tempering(self, scaling_value, input_file, output_file):
        """
        private

        Creates the scaled topology files with plumed
        """

        plumed = path.absolute_programpath(program = "plumed")

        command = [plumed, "partial_tempering", f"{scaling_value}"]

        print(f"Running Plumed, with scaling value {scaling_value}")

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

        string = "#!/bin/bash\n \n##THIS SCRIPT CREATES THE TPR FILES RUN IT BEFORE THE WORKLOADMANAGER ONE\nIT IS FOR A PLUMED PATCHED GROMACS INSTALLATION\n\n\n"

        for i in range(self.BATTERIES):
            for j in range(self.replicas):
                string = string + f"gmx grompp -maxwarn 100 -o BATTERY{i}/scaled{j}/{self.output_tpr_file} -f {self.mdp_file} -p {scaled_topologies[j]} -c {self.Protein.gro_file.rsplit('/', 1)[-1]} \n"


        write_on_files.write_file(lines = [string], file_name = self.HREM_dir + "/" + filename)

        

    def _get_hamiltonian_scaling_values(self):
        """
        Scales the hamiltonian with a geometrical progression
        scale(m) =scale^(m/(nprocs−1)) with scale = 0.2 and nprocs = self.replicas
        0 <= m <= nprocs -1

        for more information check the orac manual http://www.chim.unifi.it/orac/orac-manual.pdf
        page 123 """

        nprocs = self.replicas
        scale = 0.2

        #instantiating the list and putting the value for scale^0 = 1.0
        hamiltonian_scaling_values = [1.0]

        for m in range(1, nprocs):

            scale_m = scale ** ( m / (nprocs -1) )

            hamiltonian_scaling_values.append(scale_m)

        hamiltonian_scaling_values = tuple(hamiltonian_scaling_values)

        return hamiltonian_scaling_values


    def _change_absolute_ligand_itp_include_in_relative(self, top_file):
        """
        private

        changes the #include of the ligand itp file in the protein top from the absolute one in the relative one
        in the given topology file in order to get a meaningful #include in the self.HREM_dir directory
        that will be copied on a different Computer (the HPC cluster)
        """

        lines = read_file.read_file(file_name = top_file)

        ligand_itp_files = []

        for lgand in self.Protein.get_ligand_list():
            ligand_itp_files.append(lgand.itp_file)

        for i in range(len(lines)):

            if lines[i].strip()[:8] == "#include":

                if lines[i].split()[1].replace('"', '').strip() in ligand_itp_files:

                    itp_file = lines[i].split()[1].replace('"', '').strip()
                    itp_file = itp_file.split("/")[-1]

                    lines[i] = f'#include "{itp_file}"'

        write_on_files.write_file(lines = lines, file_name = top_file)




    def execute(self):
        """This method does not run gromacs but
        creates the input to make a REM simulation on a
        HPC cluster"""

        #crate an empty file for plumed
        write_on_files.write_file(lines = ['\n'], file_name = "empty_plumed.dat")

        #creates the REM directory that will be copied to the HPC cluster
        os.makedirs(self.HREM_dir, exist_ok=True)

        #create the BATTERY dirs and the scaled dirs
        for i in range(self.BATTERIES):

            for j in range(self.replicas):

                os.makedirs(self.HREM_dir + "/" + f"BATTERY{i}" + "/" + f"scaled{j}", exist_ok=True)

                #copy the empty file fo plumed
                shutil.copy("empty_plumed.dat", self.HREM_dir + "/" + f"BATTERY{i}" + "/" + f"scaled{j}")

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
                                        output_file = self.HREM_dir + "/" + f"{self.Protein.protein_id}_scaled_{counter}.top")

            scaled_topologies.append(f"{self.Protein.protein_id}_scaled_{counter}.top")


        #copy the needed protein files in the new directories
        for j in (self.Protein.gro_file, self.Protein.top_file):

            shutil.copy(j, self.HREM_dir)

        #copy the ligand itp files
        for lgand in self.Protein.get_ligand_list():

            shutil.copy(lgand.itp_file, self.HREM_dir)

        

        #write workload manager input for different workload managers (slurm pbs ...)
        #already in the self.HREM_dir
        self._write_workloadmanager_inputs()


        #make and copy the script that will make the tpr files in loco
        self._make_TPR_files_script(scaled_topologies = scaled_topologies)


        #create a file with important info for post processing of the HREM output
        # key = value
        important_info = [
            f"ligand_resname = {self.Protein.get_ligand_list()[0].resname}\n", #only takes the first ligand, because actually there should be ony one
            f"ligand_resnum = {self.Protein.get_ligand_list()[0].resnum}\n",
            f"ligand_itp = {self.Protein.get_ligand_list()[0].itp_file.split('/')[-1]}\n",
            f"protein_topology = {self.Protein.top_file.split('/')[-1]}\n"
        ]

        write_on_files.write_file(lines = important_info, file_name = self.HREM_dir + "/" + "important_info.dat")


        #as gmx2pdb makes some random but needed itp files some times
        #I lazily copy them all
        for file_name in os.listdir(os.getcwd()):

            if file_name[-4:] == ".itp":
                shutil.copy(file_name, self.HREM_dir)

        #change the absolute #include in relative ones as they would make no sense when the directory will be copied
        #on a different computer (the HPC cluster)
        for file_name in os.listdir(self.HREM_dir):

            if file_name[-4:] == ".top":
                self._change_absolute_ligand_itp_include_in_relative(top_file = self.HREM_dir + "/" + file_name)

        return self.Protein



class GromacsHREMOnlyLigand(GromacsHREMInput):
    """
    Makes the HREM to a system where there is only the ligand
    for any ligand in Protein._ligands

    Solvent can or cannot be present (if you want to use the output for FS-DAM there should be no solvent)
    """

    def __init__(self,
                Protein,
                only_solvent_box_gro,
                only_solvent_box_top,
                MD_program_path = 'gmx',
                kind_of_processor = 'skylake',
                number_of_cores_per_node = 64,
                use_gpu = 'auto',
                gpus_per_node = 1,
                number_of_replicas = 8):

        super().__init__(Protein = Protein,
                        MD_program_path = MD_program_path,
                        kind_of_processor = kind_of_processor,
                        number_of_cores_per_node = number_of_cores_per_node,
                        use_gpu = use_gpu,
                        gpus_per_node = gpus_per_node)

        self.only_solvent_box_gro = only_solvent_box_gro
        self.only_solvent_box_top = only_solvent_box_top

        self.BATTERIES = self._get_BATTERIES()
        #the replicas for BATTERY
        self.replicas = number_of_replicas

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


 
    def _get_BATTERIES(self):

        #for the ligand one battery is more than enough
        return 1


    def _make_TPR_files_script(self, scaled_topologies, Ligand):
        """
        private
        
        writes a bash script to create the .tpr files in loco on the HPC cluster
        it makes both the one for a plumed patched gromacs installation and for a non patched one
        """

        filename = "MAKE_TPR_FILES.sh"

        string = "#!/bin/bash\n \n##THIS SCRIPT CREATES THE TPR FILES RUN IT BEFORE THE WORKLOADMANAGER ONE\nIT IS FOR A PLUMED PATCHED GROMACS INSTALLATION\n\n\n"

        for i in range(self.BATTERIES):
            for j in range(self.replicas):
                string = string + f"gmx grompp -maxwarn 100 -o BATTERY{i}/scaled{j}/{self.output_tpr_file} -f {self.mdp_file} -p {scaled_topologies[j]} -c {Ligand.gro_file.rsplit('/', 1)[-1]} \n"


        write_on_files.write_file(lines = [string], file_name = self.HREM_dir + "/" + filename)



    def execute(self):
        """This method does not run gromacs but
        creates the input to make a REM simulation on a
        HPC cluster"""

        if self.Protein.get_ligand_list() == []:
            return self.Protein

        #crate an empty file for plumed
        write_on_files.write_file(lines = ['\n'], file_name = "empty_plumed.dat")

        Ligand = self.Protein.get_ligand_list()

        for i in range(len(Ligand)):

            #the input directory that will be created
            self.HREM_dir = f"{Ligand[i].resname}_only_ligand_HREM"

            

            self.elaborated_top_file = f"{Ligand[i].resname}_only_ligand_elaborated_topology.top"

            self.mdp_file = f"{Ligand[i].resname}_only_ligand_HREM.mdp"

            self.output_tpr_file = f"HREM.tpr"

            #creates the REM directory that will be copied to the HPC cluster
            os.makedirs(self.HREM_dir, exist_ok=True)

            #create the BATTERY dirs and the scaled dirs
            for i in range(self.BATTERIES):

                for j in range(self.replicas):

                    os.makedirs(self.HREM_dir + "/" + f"BATTERY{i}" + "/" + f"scaled{j}", exist_ok=True)

                    #copy the empty file for plumed
                    shutil.copy("empty_plumed.dat", self.HREM_dir + "/" + f"BATTERY{i}" + "/" + f"scaled{j}")

            #write the mdp file
            self._write_template_on_file()
            #and copy it in the HREM dir
            shutil.copy(self.mdp_file, self.HREM_dir)

            #create the elaborated top file for plumed
            self._interact_with_gromacs(string = [f"{self.MD_program_path}", "grompp", "-f", f"{self.mdp_file}", "-c", f"{Ligand[i].gro_file}", "-p", f"{Ligand[i].top_file}", "-maxwarn", "100", "-pp", f"{self.elaborated_top_file}"])

            #finds the hot residues and edits the elaborated top file accordingly
            self._edit_top_file(top_file = self.elaborated_top_file)


            #use plumed to make the scaled topologies

            hamiltonian_scaling_values = self._get_hamiltonian_scaling_values()
            scaled_topologies = []
            for counter, value in enumerate(hamiltonian_scaling_values):
                #they are already saved in the HREM dir
                self._plumed_partial_tempering(scaling_value = value,
                                            input_file = self.elaborated_top_file,
                                            output_file = self.HREM_dir + "/" + f"{Ligand[i].resname}_scaled_{counter}.top")

                scaled_topologies.append(f"{Ligand[i].resname}_scaled_{counter}.top")


            #copy the needed files in the new directories
            for j in (Ligand[i].gro_file, Ligand[i].top_file, Ligand[i].itp_file, self.only_solvent_box_gro, self.only_solvent_box_top):

                shutil.copy(j, self.HREM_dir)

            
            #write workload manager input for different workload managers (slurm pbs ...)
            #already in the self.HREM_dir
            self._write_workloadmanager_inputs()
            

            #make and copy the script that will make the tpr files in loco
            self._make_TPR_files_script(scaled_topologies = scaled_topologies, Ligand = Ligand[i])


            #create a file with important info for post processing of the HREM output
            # key = value
            important_info = [
                f"ligand_resname = {Ligand[i].resname}\n",
                f"ligand_itp = {Ligand[i].itp_file.split('/')[-1]}\n",
                f"ligand_top = {Ligand[i].top_file.split('/')[-1]}\n",
                f"only_solvent_gro = {self.only_solvent_box_gro.split('/')[-1]}\n",
                f"only_solvent_top = {self.only_solvent_box_top.split('/')[-1]}\n"
            ]

            write_on_files.write_file(lines = important_info, file_name = self.HREM_dir + "/" + "important_info.dat")


            #as gmx2pdb makes some random but needed itp files some times
            #I lazily copy them all
            for file_name in os.listdir(os.getcwd()):

                if file_name[-4:] == ".itp":
                    shutil.copy(file_name, self.HREM_dir)

            #change the absolute #include in relative ones as they would make no sense when the directory will be copied
            #on a different computer (the HPC cluster)
            for file_name in os.listdir(self.HREM_dir):

                if file_name[-4:] == ".top":
                    self._change_absolute_ligand_itp_include_in_relative(top_file = self.HREM_dir + "/" + file_name)



        return self.Protein





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
            
            string = string + f"mpirun -np {cpus_per_task * self.replicas} gmx_mpi mdrun {gpu_options(use_gpu = self.gpu)} -v -plumed empty_plumed.dat -replex 100 -hrex -dlb no {multidir_string(i)} -s {self.tpr_file} & \n"

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


    