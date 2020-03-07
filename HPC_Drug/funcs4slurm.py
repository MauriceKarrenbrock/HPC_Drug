"""Contains the functions to make a slurm input (workload manager)"""
from HPC_Drug import pipeline_functions


import math

class SlurmInput(object):
    """This clas supports other classes writing the
    slurm input to make the other classes' inputs run
    on HPC clusters
    """

    def __init__(self,
                MD_input_file = 'md.in',
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
                MD_program_path = None):


        #if it is a string I transform it in a one object list
        #to iterate on it if necessary
        if type(MD_input_file) == str:
            MD_input_file = [MD_input_file]        
        self.MD_input_file = MD_input_file


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

        string = ""
        for i in pipeline_functions.get_iterable(self.MD_input_file):
            string = string + f"mpirun 	  {self.MD_program_path} < {i}\n"

        return string


    def write_gromacs_mpirun_string(self):
        """private"""

        string = ""

        for i in pipeline_functions.get_iterable(self.MD_input_file):
            string = string + f"mpirun -np {self.cpus_per_task * 8} gmx_mpi mdrun_d -v -plumed empty_plumed.dat -multi 8 -replex 100 -hrex -dlb no -s {i} & \n"

        string = string + "wait"

        return string

    def write(self):
        """Writes the chosen template on a file called as 
        self.slurm_input_file """

        with open(self.slurm_input_file, 'w') as w:

            for line in self.template:

                w.write(f"{line.strip()}\n")

        return self.slurm_input_file
