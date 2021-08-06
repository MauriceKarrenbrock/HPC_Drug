######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
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


parser = argparse.ArgumentParser(
    description="This script creates the input for FS-DAM both for the boded and unbonded system. "
    "For Gromacs\n You must run this script in the HREM root directory",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--program-path',
    action = "store",
    default = "gmx",
    type=str,
    help = "The absolute path to the chosen program executable")

parser.add_argument('--hrem-type',
    action="store",
    default="protein-ligand",
    type=str,
    choices = ["protein-ligand", "only-ligand"],
    help = "If the HREM dir contains the output of the ligand alone (unbound) or the protein-ligand system (bound) (default)")

parser.add_argument('--vdw-timestep-ps',
    action="store",
    default=None,
    type=float,
    help = "The timestep for the vdw transformation, the default will work most of the times")

parser.add_argument('--q-timestep-ps',
    action="store",
    default=None,
    type=float,
    help = "The timestep for the q transformation, the default will work most of the times")

parser.add_argument('--vdw-number-of-steps',
    action="store",
    default=None,
    type=int,
    help = "The number of steps for the q transformation, the default will work most of the times")

parser.add_argument('--q-number-of-steps', action="store", default=None, help = "The number of steps for the q transformation, the default will work most of the times")

parser.add_argument('--pbc-atoms',
    action="store",
    default=None,
    type=str,
    help = "The pbc atoms (see gromacs documentation), COMMA SEPARATED LIST, are only needed when making a protein-ligand transformation, the program will try to get them autonomously but it might be unable to do it in complex situations")

parser.add_argument('--number-of-frames-to-use',
    action="store",
    default=200,
    type=int,
    help = "The number of alchemical transformations to do, a good default is 200 protein ligand and 400 ligand only, default=200")

parser.add_argument('--constrains',
    action="store",
    default=None,
    type=str,
    choices = ['none', 'h-bonds', 'all-bonds'],
    help = "how to constraints vibrations, same syntax as gromacs mdp file default=all-bonds")

parser.add_argument('--reference-frame',
    action="store",
    default=None,
    type=str,
    help = "A PDB or GRO file to use as reference to get the nearest neighboring residues. "
    "and if needed the number of atoms in the binding pocket "
    "if not givent the first frame in the BATTERY0 trajectory will be used")

parser.add_argument('--atoms-in-pocket',
    action="store",
    default=None,
    type=int,
    help = "The number of atoms that must be in the binding pocket uìin order to consider the ligand "
    "in the pocket. If not given the number of atoms in the --reference-frame will be used ")


parsed_input = parser.parse_args()

parsed_input.program_path = path.absolute_programpath(parsed_input.program_path)

if parsed_input.hrem_type == "protein-ligand":

    creation=False

elif parsed_input.hrem_type == "only-ligand":

    creation=True

# From comma separated list to list of int
if parsed_input.pbc_atoms is not None:
    parsed_input.pbc_atoms = parsed_input.pbc_atoms.split(',')
    for i, atom in enumerate(parsed_input.pbc_atoms):
        parsed_input.pbc_atoms[i] = int(atom)

fsdam_obj = fsdam.FSDAMInputPreprocessing(gromacs_path = parsed_input.program_path,
                vdw_timestep_ps=parsed_input.vdw_timestep_ps,
                q_timestep_ps=parsed_input.q_timestep_ps,
                vdw_number_of_steps=parsed_input.vdw_number_of_steps,
                q_number_of_steps=parsed_input.q_number_of_steps,
                creation=creation,
                pbc_atoms=parsed_input.pbc_atoms,
                number_of_frames_to_use=parsed_input.number_of_frames_to_use,
                constrains=parsed_input.constrains,
                reference_frame=parsed_input.reference_frame,
                atoms_in_pocket=parsed_input.atoms_in_pocket)

fsdam_obj.execute()
