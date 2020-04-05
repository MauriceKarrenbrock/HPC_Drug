######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains the functions to make a slurm input (workload manager)
"""

from HPC_Drug import pipeline_functions
from HPC_Drug.auxiliary_functions import get_iterable

import collections
import math

class SlurmInput(object):
    """This clas supports other classes writing the
    slurm input to make the other classes' inputs run
    on HPC clusters
    """

    def __init__(self,
                MD_input_file = 'md.in',
                MD_input_directories = None,
                slurm_input_file = 'HPC_Drug.slr',
                MD_program = 'gromacs',
                MD_calculation_type = 'rem',
                number_of_cores_per_node = 64,
                max_time = "24:00:00",
                ntasks = 64,
                cpus_per_task = 8,
                std_out = 'STD_out.out',
                std_err = 'STD_err.err',
                partition_name = None,
                account_name = None,
                MD_program_path = None,
                use_gpu = 'auto'):

        
        self.MD_input_file = MD_input_file


        #it should be a list, it can be a nested list with 2 indexes
        #the first will be the mpirun index, the second will be the one representing the
        #directories to give to gmx -multidir (works only for gromacs not orac)
        #example:
        #self.MD_input_directories = [['a', 'b'], ['c', 'd']]
        #/usr/bin/bash
        # ....
        #mpirun -multidir a b &
        #mpirun -multidir c d &
        #wait
        self.MD_input_directories = MD_input_directories
        
        if self.MD_input_directories == None:
            pass
        elif type(self.MD_input_directories) == str:
            #if it is a string I adapt the file name to it
            #it is actually a quite bad way of using it
            #it would be better to directly give the right MD_input_file and keep it None
            self.MD_input_file = self.MD_input_directories.strip().rstrip('/') + '/' + self.MD_input_file.strip().lstrip('/')

            self.MD_input_directories = None

        #if I was given a monodimensional list I transform it in a bidimensional one
        elif isinstance(self.MD_input_directories[0], collections.Iterable) == False:

            self.MD_input_directories = [self.MD_input_directories]




        self.slurm_input_file = slurm_input_file
        self.MD_program = MD_program.lower().strip()
        self.MD_calculation_type = MD_calculation_type.lower().strip()

        self.number_of_cores_per_node = number_of_cores_per_node

        #the maximum amount of wall time alowed hh:mm:ss
        if type(max_time) == int:
            max_time = str(max_time)
        elif len(max_time.split(":")) > 3 or len(max_time.split(":")) == 0:
            raise ValueError(f"{max_time} is not a valid max_time value, must be hh:mm:ss")
        while len(max_time.split(":")) < 2:
            max_time = max_time + ":00"
        self.max_time = max_time


        self.ntasks = ntasks

        self.cpus_per_task = cpus_per_task
        self.std_out = std_out
        self.std_err = std_err


        #many HPC clusters have different partitions
        self.partition_name = partition_name

        self.account_name = account_name
        
        self.MD_program_path = MD_program_path
        if self.MD_program_path == None and self.MD_program == "gromacs":
            self.MD_program_path = "gmx"

        elif self.MD_program_path == None and self.MD_program == "orac":
            self.MD_program_path = "~/ORAC/trunk/src/INTEL-FFTW-OMP-MPI/orac"

        #gromacs has various options to use gpu
        #auto (default) that will use all the available ones automaticly
        #cpu uses no GPU even if available
        #gpu forces the use of GPU (but in case you want to use a gpu auto would be safer and more robust)
        self.use_gpu = use_gpu.lower().strip()
        if self.use_gpu not in ('auto', 'cpu', 'gpu'):
            raise ValueError(f"{self.use_gpu} is not a valid gpu option, valid options are auto cpu gpu")

        self.template = self.get_template()

    def get_template(self):
        """private function used by the __init__ function
        to get the right template"""

        REM_SBATCH_template = [
            "#!/bin/bash",
            "\n##CAUTION!!!!",
            "##THIS IS A VERY GENERIC SLURM INPUT FILE",
            "##CONTACT YOUR SISTEM ADMINISTRATOR TO KNOW WHICH EXTRA LINES YOU MAY NEED TO INSERT HERE",

            self.write_account_string(),
            
            f"#SBATCH --time={self.max_time}",

            f"#SBATCH --nodes={math.ceil((self.ntasks) / (math.floor((self.number_of_cores_per_node) / (self.cpus_per_task))))} # ntasks / ntasks per node",
            
            f"#SBATCH --ntasks-per-node={math.floor(self.number_of_cores_per_node / self.cpus_per_task)} # math.floor(cpu_per_node / cpus_per_task)",
            
            f"#SBATCH --ntasks={self.ntasks} ## BATTERIES * 8 (given as input)",
            
            f"#SBATCH --cpus-per-task={self.cpus_per_task} ##usually 8",

            self.write_partition_string(),
            
            f"#SBATCH --output={self.std_out}",
            
            f"#SBATCH --error={self.std_err}",

            " "]

        orac_REM_mpirun_template =[
            self.write_orac_mpirun_string()
        ]

        gromacs_REM_mpirun_template = [
            self.write_gromacs_mpirun_string()
        ]

        if self.MD_program == 'orac' and self.MD_calculation_type == 'rem':

            orac_REM_template = REM_SBATCH_template + orac_REM_mpirun_template 

            return orac_REM_template

        elif self.MD_program == 'gromacs' and self.MD_calculation_type == 'rem':

            gromacs_REM_template = REM_SBATCH_template + gromacs_REM_mpirun_template 

            return gromacs_REM_template

    def write_partition_string(self):
        """private"""

        if self.partition_name == None:
            return ''

        else:
            return f"#SBATCH --partition={self.partition_name}"

    def write_account_string(self):
        """private"""

        if self.account_name == None:
            return " "

        else:
            return f"#SBATCH -A  {self.account_name}"

    def write_orac_mpirun_string(self):
        """private"""

        if type(self.MD_input_file) == str:

            string = f"mpirun 	  {self.MD_program_path} < {self.MD_input_file}\n"

        else :

            string = ""
            
            for i in get_iterable.get_iterable(self.MD_input_file):
                string = string + f"mpirun 	  {self.MD_program_path} < {i}\n"

        return string


    def write_gromacs_mpirun_string(self):
        """private"""

        def gpu_options(use_gpu):
            if use_gpu == 'auto':
                return ""

            else:
                string = f" -nb {use_gpu} -pme {use_gpu} -bonded {use_gpu} -update {use_gpu} -pmefft {use_gpu}"
                
                if use_gpu == 'gpu':
                    string = string + " -npme 1"

                return string

        def multidir_string(dirs):
            
            string = " -multidir "

            for i in dirs:
                string = string + f" {i} "

            return string


            

        string = ""

        if self.MD_input_directories == None:

            string = f"mpirun -np {self.cpus_per_task * 8} gmx_mpi mdrun {gpu_options(use_gpu = self.use_gpu)} -v -plumed empty_plumed.dat -replex 100 -hrex -dlb no -s {self.MD_input_file} \n"

        else:
            for dirs in self.MD_input_directories:
                string = string + f"mpirun -np {self.cpus_per_task * 8} gmx_mpi mdrun {gpu_options(use_gpu = self.use_gpu)} -v -plumed empty_plumed.dat -replex 100 -hrex -dlb no {multidir_string(dirs = dirs)} -s {self.MD_input_file} & \n"

            string = string + "wait"

        return string

    def write(self):
        """Writes the chosen template on a file called as 
        self.slurm_input_file """

        with open(self.slurm_input_file, 'w') as w:

            for line in self.template:

                w.write(f"{line.strip()}\n")

        return self.slurm_input_file
