# -*- coding: utf-8 -*-
# pylint: disable=line-too-long
# pylint: disable=duplicate-code
#############################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""script to do post processing of a FSDAM or a vDSSB

use --help for usage info
"""

import argparse
import shutil
from pathlib import Path

from PythonFSDAM.pipelines import superclasses

from HPC_Drug.MD.gromacs.volume_correction import volume_correction as gromacs_volume_correction

parser = argparse.ArgumentParser(
    description='This script will post process everything needed after '
    'FSDAM (vDSSB) use --help for usage info '
    'Parallelism by default will use all the cores evailable, '
    'use OMP_NUM_THREADS environment variable to limit '
    'the number of used CPUs',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--md-program',
                    action='store',
                    type=str,
                    default='gromacs',
                    choices=['gromacs'],
                    help='The MD program to use')

parser.add_argument(
    '--unbound-dir',
    action='store',
    type=str,
    help=
    'the directory in which the various RESTART and Extra_RESTART* directories are for '
    'the ligand system, the program will use them and detect them all')

parser.add_argument(
    '--bound-dir',
    action='store',
    type=str,
    help=
    'the directory in which the various RESTART and Extra_RESTART* directories are for '
    'the protein-ligand system, the program will use them and detect them all')

parser.add_argument('--temperature',
                    action='store',
                    type=float,
                    default=298.15,
                    help='temperature in Kelvin (K)')

parsed_input = parser.parse_args()


bound_dir = Path(parsed_input.bound_dir)
unbound_dir = Path(parsed_input.unbound_dir)

if parsed_input.md_program == 'gromacs':

    # Bound files
    bound_files = [[], []]

    bound_files[0] = bound_dir.glob('RESTART/q*.xvg')
    bound_files[0] = [i.resolve() for i in bound_files[0] if 'pull' not in str(i)]

    bound_files[1] = bound_dir.glob('RESTART/vdw*.xvg')
    bound_files[1] = [i.resolve() for i in bound_files[1] if 'pull' not in str(i)]

    i = 0
    while (bound_dir / f'Extra_RESTART{i}').exists():
        tmp_files = bound_dir.glob(f'Extra_RESTART{i}/q*.xvg')
        bound_files[0] += [i.resolve() for i in tmp_files if 'pull' not in str(i)]

        tmp_files = bound_dir.glob('RESTART/vdw*.xvg')
        bound_files[1] += [i.resolve() for i in tmp_files if 'pull' not in str(i)]

        i += 1

    # Unbound files
    unbound_files = [[], []]

    unbound_files[0] = unbound_dir.glob('RESTART/vdw*.xvg')
    unbound_files[0] = [i.resolve() for i in unbound_files[0] if 'pull' not in str(i)]

    unbound_files[1] = unbound_dir.glob('RESTART/q*.xvg')
    unbound_files[1] = [i.resolve() for i in unbound_files[1] if 'pull' not in str(i)]

    i = 0
    while (unbound_dir / f'Extra_RESTART{i}').exists():
        tmp_files = unbound_dir.glob(f'Extra_RESTART{i}/vdw*.xvg')
        unbound_files[0] += [i.resolve() for i in tmp_files if 'pull' not in str(i)]

        tmp_files = unbound_dir.glob('RESTART/q*.xvg')
        unbound_files[1] += [i.resolve() for i in tmp_files if 'pull' not in str(i)]

        i += 1

    #################################
    # Volume free energy correction
    volume_correction = gromacs_volume_correction(directory=bound_dir, temperature=parsed_input.temperature)

    with open('volume_correction.dat', 'w') as f:
        f.write('# Volume correction in Kcal/mol, must be added to the dissociation '
        'free energy or subtracted from the binding free energy\n'
        f'{volume_correction:.18e}\n')

    del volume_correction

    print('created volume_correction.dat with the volume correction')
    ##################################


if len(bound_files) == 1:

    bound_files = bound_files[0]

print('Calculating Jarzinski free energy, will take some time')
obj = superclasses.JarzynskiVDSSBPostProcessingPipeline(
    bound_files,
    unbound_files,
    temperature=parsed_input.temperature,
    md_program=parsed_input.md_program)

free_energy, std = obj.execute()

print(f'Jarzynski free energy {free_energy:.18e} Kcal/mol\nCI95 {std*1.96:.18e}')

print('Calculating Gaussian Mixtures (EM) free energy, will take some time')
obj = superclasses.GaussianMixturesVDSSBPostProcessingPipeline(
    bound_files,
    unbound_files,
    temperature=parsed_input.temperature,
    md_program=parsed_input.md_program)

free_energy, std = obj.execute()

print(f'Gaussian mixtures (EM) free energy {free_energy:.18e} Kcal/mol\nCI95 {std*1.96:.18e}')

