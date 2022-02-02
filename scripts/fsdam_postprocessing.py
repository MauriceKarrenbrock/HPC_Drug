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
from pathlib import Path
import shutil

from PythonFSDAM.pipelines import superclasses
from PythonFSDAM import integrate_works

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

parser.add_argument('--kind-of-process',
                    action='store',
                    type=str,
                    default='vdssb',
                    choices=['fsdam-unbound-creation', 'fsdam-unbound-annihilation', 'vdssb'],
                    help='if you are doing a standard FSDAM or a vDSSB, for the bound system (protein + ligand) '
                    'it will always be taken for granted that the ligand was annihilated, if you are doing vDSSB the ligand MUST '
                    'be created in a box of water, if you are doing standard FSDAM you have to specify what you are doing '
                    'with the unbound system')

parsed_input = parser.parse_args()


if parsed_input.kind_of_process == 'fsdam-unbound-annihilation':
    ligand_creation = False
else:
    ligand_creation = True


bound_dir = Path(parsed_input.bound_dir)
unbound_dir = Path(parsed_input.unbound_dir)

if parsed_input.md_program == 'gromacs':

    # Bound files
    bound_files = [[], []]

    bound_files[0] = bound_dir.glob('RESTART/q*.xvg')
    bound_files[0] = sorted([i.resolve() for i in bound_files[0] if 'pull' not in str(i)])

    bound_files[1] = bound_dir.glob('RESTART/vdw*.xvg')
    bound_files[1] = sorted([i.resolve() for i in bound_files[1] if 'pull' not in str(i)])

    i = 0
    while (bound_dir / f'Extra_RESTART{i}').exists():
        tmp_files = bound_dir.glob(f'Extra_RESTART{i}/q*.xvg')
        bound_files[0] += sorted([i.resolve() for i in tmp_files if 'pull' not in str(i)])

        tmp_files = bound_dir.glob('RESTART/vdw*.xvg')
        bound_files[1] += sorted([i.resolve() for i in tmp_files if 'pull' not in str(i)])

        i += 1

    # Unbound files
    unbound_files = [[], []]

    unbound_files[0] = unbound_dir.glob('RESTART/vdw*.xvg')
    unbound_files[0] = sorted([i.resolve() for i in unbound_files[0] if 'pull' not in str(i)])

    unbound_files[1] = unbound_dir.glob('RESTART/q*.xvg')
    unbound_files[1] = sorted([i.resolve() for i in unbound_files[1] if 'pull' not in str(i)])

    i = 0
    while (unbound_dir / f'Extra_RESTART{i}').exists():
        tmp_files = unbound_dir.glob(f'Extra_RESTART{i}/vdw*.xvg')
        unbound_files[0] += sorted([i.resolve() for i in tmp_files if 'pull' not in str(i)])

        tmp_files = unbound_dir.glob('RESTART/q*.xvg')
        unbound_files[1] += sorted([i.resolve() for i in tmp_files if 'pull' not in str(i)])

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


################################################
# Do vDSSB
if parsed_input.kind_of_process == 'vdssb':

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

##########################################################
# Do normal FSDAM
else:

    # Jarzynski
    print('Calculating Jarzinski free energy, will take some time')
    unbound_obj = superclasses.JarzynskiPostProcessingAlchemicalLeg(
        unbound_files,
        temperature=parsed_input.temperature,
        md_program=parsed_input.md_program,
        creation=ligand_creation)

    unbound_free_energy, unbound_std = unbound_obj.execute()

    shutil.move(f'{str(unbound_obj)}_free_energy.dat',
                'unbound_' + f'{str(unbound_obj)}_free_energy.dat')

    shutil.move('work_values.dat',
                'unbound_work_values.dat')

    print(f'Jarzynski unbound free energy {unbound_free_energy}\n' f'CI95 {1.96*(unbound_std)}')

    bound_obj = superclasses.JarzynskiPostProcessingAlchemicalLeg(
        bound_files,
        temperature=parsed_input.temperature,
        md_program=parsed_input.md_program,
        creation=False)

    bound_free_energy, bound_std = bound_obj.execute()

    shutil.move(f'{str(bound_obj)}_free_energy.dat',
                'bound_' + f'{str(bound_obj)}_free_energy.dat')

    shutil.move('work_values.dat',
                'bound_work_values.dat')

    print(f'Jarzynski bound free energy {bound_free_energy}\n' f'CI95 {1.96*(bound_std)}')

    if ligand_creation:
        total_free_energy = bound_free_energy + unbound_free_energy
    else:
        total_free_energy = bound_free_energy - unbound_free_energy

    total_std = bound_std + unbound_std

    with open(f'{str(bound_obj)}_total_free_energy.dat', 'w') as f:

        f.write(
            '#total unbinding free energy (no volume correction done) in Kcal mol\n'
            '#Dg   CI95%  STD\n'
            f'{total_free_energy:.18e} {1.96*(total_std):.18e} {total_std:.18e}\n'
        )

    print(f'Jarzynski total free energy {total_free_energy}\n' f'CI95 {1.96*(total_std)}')


    # EM gaussian mixtures
    print('Calculating Gaussian mixtures (EM) free energy, will take some time')
    unbound_obj = superclasses.GaussianMixturesPostProcessingAlchemicalLeg(
        unbound_files,
        temperature=parsed_input.temperature,
        md_program=parsed_input.md_program,
        creation=ligand_creation)

    unbound_free_energy, unbound_std = unbound_obj.execute()

    shutil.move(f'{str(unbound_obj)}_free_energy.dat',
                'unbound_' + f'{str(unbound_obj)}_free_energy.dat')

    shutil.move('work_values.dat',
                'unbound_work_values.dat')

    print(f'Gaussian mixtures (EM) unbound free energy {unbound_free_energy}\n' f'CI95 {1.96*(unbound_std)}')

    bound_obj = superclasses.GaussianMixturesPostProcessingAlchemicalLeg(
        bound_files,
        temperature=parsed_input.temperature,
        md_program=parsed_input.md_program,
        creation=False)

    bound_free_energy, bound_std = bound_obj.execute()

    shutil.move(f'{str(bound_obj)}_free_energy.dat',
                'bound_' + f'{str(bound_obj)}_free_energy.dat')

    shutil.move('work_values.dat',
                'bound_work_values.dat')

    print(f'Gaussian mixtures (EM) bound free energy {bound_free_energy}\n' f'CI95 {1.96*(bound_std)}')

    if ligand_creation:
        total_free_energy = bound_free_energy + unbound_free_energy
    else:
        total_free_energy = bound_free_energy - unbound_free_energy

    total_std = bound_std + unbound_std

    with open(f'{str(bound_obj)}_total_free_energy.dat', 'w') as f:

        f.write(
            '#total unbinding free energy (no volume correction done) in Kcal mol\n'
            '#Dg   CI95%  STD\n'
            f'{total_free_energy:.18e} {1.96*(total_std):.18e} {total_std:.18e}\n'
        )

    print(f'Gaussian mixtures (EM) total free energy {total_free_energy}\n' f'CI95 {1.96*(total_std)}')

##########################################################

print('Creating csv files with work vs lambda (useful for plotting)')
# Bound
bound_csv = integrate_works.make_work_vs_lambda_csv(work_files=bound_files,
                            md_program=parsed_input.md_program,
                            creation=False,
                            num_runs=len(bound_files))

shutil.move(bound_csv,
            f'bound_{bound_csv}')

# Unbound
unbound_csv = integrate_works.make_work_vs_lambda_csv(work_files=unbound_files,
                            md_program=parsed_input.md_program,
                            creation=ligand_creation,
                            num_runs=len(bound_files))

shutil.move(unbound_csv,
            f'unbound_{bound_csv}')

print(f'Created {bound_csv} and {unbound_csv}')
