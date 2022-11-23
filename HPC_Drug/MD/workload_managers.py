######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
Contains the functions to create the header to ask for resouces allocation for many workload managesr (es slurm, pdb)
"""

class WorkloadManagerHeader(object):
    """
    The super class for all the workload manager headers
    """

    def __init__(self,
                nodes = 1,
                tasks = 1,
                tasks_per_node = 1,
                cpus_per_task = 1,
                GPUs = None,
                wall_time = "24:00:00",
                output = "output.out",
                error = "error.err",
                account_name = None,
                partition_name = None):
        """
        nodes :: integer, default 1, number of nodes to ask 
        tasks :: integer, total number of processes, default 1 
        tasks_per_node :: integer, number of processes per node (NOT the number of cores), default 1 
        cpus_per_task :: integer, number of cores per task, default 1
        GPUs :: integer, total number of GPU, default None, omitted 
        wall_time :: string, the maximum job time in hh:mm:ss, default 24:00:00 
        output :: string, default output.out, the file where to print the standard output
        error :: string, default error.err, the file where to print the standard error
        account_name :: string, your account name, default None (omitted)
        partition_name :: string, the partition that you are going to use, default None (omitted)
        """

        self.nodes = nodes
        self.tasks = tasks
        self.tasks_per_node = tasks_per_node
        self.cpus_per_task = cpus_per_task
        self.GPUs = GPUs
        self.wall_time = wall_time
        self.output = output
        self.error = error
        self.account_name = account_name
        self.partition_name = partition_name


    def _get_header(self):
        return []


    def execute(self):

        return self._get_header()



class SlurmHeader(WorkloadManagerHeader):
    """The execute() method returns a list of strings (with no newline characters, you shall add them later es '\n'.join()) that creates the #SBATCH header
    """

    
    def _get_header(self):

        slurm_header = [
            "#!/bin/bash",
            "\n##CAUTION!!!!",
            "##THIS IS A VERY GENERIC SLURM INPUT FILE",
            "##CONTACT YOUR SISTEM ADMINISTRATOR TO KNOW WHICH EXTRA LINES YOU MAY NEED TO INSERT HERE",

            self._write_account_string(),
            
            f"#SBATCH --time={self.wall_time}",

            f"#SBATCH --nodes={self.nodes}",

            f"#SBATCH --ntasks={self.tasks}",
            
            f"#SBATCH --ntasks-per-node={self.tasks_per_node}",
            
            f"#SBATCH --cpus-per-task={self.cpus_per_task}",

            self._write_gpu_string(),

            self._write_partition_string(),
            
            f"#SBATCH --output={self.output}",
            
            f"#SBATCH --error={self.error}"

            " ",
        ]

        return slurm_header

    def _write_account_string(self):
        """private"""

        if self.account_name is None:
            return ""

        else:
            return f"#SBATCH -A  {self.account_name}"

    def _write_gpu_string(self):

        if self.GPUs is None:
            return ""

        return f"#SBATCH --gres=gpu:{self.GPUs}"

    def _write_partition_string(self):
        """private"""

        if self.partition_name is None:
            return ''

        else:
            return f"#SBATCH --partition={self.partition_name}"



class PBSHeader(WorkloadManagerHeader):
    """The execute() method returns a list of strings (with no newline characters, you shall add them later es '\n'.join()) that creates the #SBATCH header
    """

    
    def _get_header(self):

        pbs_header = [
            "#!/bin/bash",
            "\n##CAUTION!!!!",
            "##THIS IS A VERY GENERIC PBS INPUT FILE",
            "##CONTACT YOUR SISTEM ADMINISTRATOR TO KNOW WHICH EXTRA LINES YOU MAY NEED TO INSERT HERE",

            f"#PBS -l walltime={self.wall_time}",

            f"#PBS -l select={self.nodes}:ncpus={self.tasks * self.cpus_per_task}:mpiprocs={self.tasks}{self._write_gpu_string()}",

            f"#PBS -o {self.output}",

            f"#PBS -e {self.error}",

            "cd $PBS_O_WORKDIR",
            " "
        ]

        return pbs_header


    def _write_gpu_string(self):

        if self.GPUs is None:
            return ""

        return f":ngpus={self.GPUs}"


