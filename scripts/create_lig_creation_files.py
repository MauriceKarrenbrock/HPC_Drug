######################################################################################
# Copyright (c) 2020-2023 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

import argparse
import mdtraj
from pathlib import Path
import HPC_Drug.MD.gromacs.fsdam as _fsdam

parser = argparse.ArgumentParser(
    description='This script will create the input files '
    'needed to create a ligand in a protein pocket '
    'it works by superimposing frames from the replica exchange of the ligand '
    'in vacuum on the binding pocket of the apo protein by using frames of the '
    'replica exchange of the holo system as a reference\n'
    'All molecules must be whole\n'
    'WARNING: This script is a bit naive and sometimes the ligand does not get put in the binding pocket '
    '(I do not know why), therefore check the output files carefully',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--output-dir',
                    action='store',
                    type=str,
                    default='.',
                    help='The directory in which the output files should be saved, '
                    'default is the current working directory')

parser.add_argument('--apo-files',
                    action='store',
                    type=str,
                    help='A comma separated list of files of the apo replica exchange run(s) '
                    'any format supported by mdtraj is accepted (xtc, trr, pdb, gro, ...)')

parser.add_argument('--apo-files-top',
                    action='store',
                    type=str,
                    help='A file of the apo system to use as topology for mdtraj (like pdb or gro files)')

parser.add_argument('--holo-files',
                    action='store',
                    type=str,
                    help='A comma separated list of files of the holo replica exchange run(s) '
                    'any format supported by mdtraj is accepted (xtc, trr, pdb, gro, ...)')

parser.add_argument('--holo-files-top',
                    action='store',
                    type=str,
                    help='A file of the holo system to use as topology for mdtraj (like pdb or gro files)')

parser.add_argument('--ligand-in-vacuum-files',
                    action='store',
                    type=str,
                    help='A comma separated list of files of the ligand in vacuum '
                    'replica exchange run(s) '
                    'any format supported by mdtraj is accepted (xtc, trr, pdb, gro, ...)')

parser.add_argument('--ligand-in-vacuum-files-top',
                    action='store',
                    type=str,
                    help='A file of the ligand in vacuum system to use as topology for '
                    'mdtraj (like pdb or gro files)')

parser.add_argument('--lig-selection-string',
                    action='store',
                    type=str,
                    help='A mdtraj selection string that allows to select the ligand')


parsed_input = parser.parse_args()

output_dir = Path(parsed_input.output_dir).resolve()
output_dir.mkdir(parents=True, exist_ok=True)


apo_files = mdtraj.load(parsed_input.apo_files.strip().split(","), top=parsed_input.apo_files_top)

holo_files = mdtraj.load(parsed_input.holo_files.strip().split(","), top=parsed_input.holo_files_top)

ligand_in_vacuum_files = mdtraj.load(parsed_input.ligand_in_vacuum_files.strip().split(","),
                                    top=parsed_input.ligand_in_vacuum_files_top)


_fsdam.create_complex_lig_creation_frames(apo_files=parsed_input.apo_files.strip().split(","),
                                    holo_files=parsed_input.holo_files.strip().split(","),
                                    ligand_in_vacuum_files=parsed_input.ligand_in_vacuum_files.strip().split(","),
                                    apo_files_top=parsed_input.apo_files_top, 
                                    holo_files_top=parsed_input.holo_files_top,
                                    ligand_in_vacuum_files_top=parsed_input.ligand_in_vacuum_files_top,
                                    lig_selection_string=parsed_input.lig_selection_string,
                                    output_dir=parsed_input.output_dir,
                                    verbose=True)
