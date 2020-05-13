######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

import os
import argparse

from HPC_Drug.auxiliary_functions import path
from HPC_Drug.MD.gromacs import fsdam


parser = argparse.ArgumentParser(description="This script creates the input for FS-DAM both for the boded and unbonded system. For Gromacs or Orac (not implemented yet)")

parser.add_argument('--program', action="store", default="gromacs", choices = ["gromacs", "orac"], help = "The program that shall be used")
parser.add_argument('--program-path', action = "store", default = "gmx", help = "The absolute path to the chosen program executable")
parser.add_argument('--hrem-dir', action="store", default=os.getcwd(), help = "The root directory of the HREM output, default current working directory")

parser.add_argument('--hrem-type', action="store", default="protein-ligand", choices = ["protein-ligand", "only-ligand"], help = "If the HREM dir contains the output of the ligand alone (unbound) or the protein-ligand system (bound) (default)")

parser.add_argument('--use-gpu', action="store_true", default=False, help = "If omitted the program will take for granted that you do not want to use GPUs (or that your HPC cluster of choice doesn't have them)")
parser.add_argument('--gpu-per-node', action="store", default=1, type = int, help = "The number of GPUs per node, default 1, will be ignored if not --use-gpu")
parser.add_argument('--cpu-per-node', action="store", default=64, type = int, help = "The number of CPUs per node, default 64")

parsed_input = parser.parse_args()

parsed_input.program_path = path.absolute_programpath(parsed_input.program_path)

parsed_input.hrem_dir = path.absolute_filepath(parsed_input.hrem_dir)



if parsed_input.program == "gromacs":

    if parsed_input.hrem_type == "protein-ligand":
    
        fsdam_obj = fsdam.FSDAMInputProteinLigand(HREM_dir = parsed_input.hrem_dir,
                                                gromacs_path = parsed_input.program_path,
                                                cpus_per_node = parsed_input.cpu_per_node,
                                                use_GPU = parsed_input.use_gpu,
                                                GPU_per_node = parsed_input.use_gpu)

    elif parsed_input.hrem_type == "only-ligand":

        fsdam_obj = fsdam.FSDAMInputOnlyLigand(HREM_dir = parsed_input.hrem_dir,
                                            gromacs_path = parsed_input.program_path,
                                            cpus_per_node = parsed_input.cpu_per_node,
                                            use_GPU = parsed_input.use_gpu,
                                            GPU_per_node = parsed_input.use_gpu)

    fsdam_obj.execute()

elif parsed_input.program == "orac":
    raise NotImplementedError("Sorry, but will be implemented soon")


