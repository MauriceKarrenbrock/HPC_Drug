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
import warnings
import traceback

import mdtraj

from PythonFSDAM.pipelines import superclasses
from PythonFSDAM import integrate_works
from PythonFSDAM import free_energy_charge_correction as _charge_correction

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
                    choices=['fsdam', 'vdssb'],
                    help='if you are doing a standard FSDAM or a vDSSB, '
                    'remember that if you are doing vDSSB if you are annihilating the ligand in the bound system you must be creatighn it in the unbound '
                    'and viceversa (use the --bound-creation and --unbound-creation flags)')

parser.add_argument('--bound-creation',
                    action='store_true',
                    help='If in the bound transformation you are creating the ligand in the protein-ligand system use this flag, '
                    'remember that if you are doing vDSSB if you are annihilating the ligand in the bound system you must be creatighn it in the unbound '
                    'and viceversa')

parser.add_argument('--unbound-creation',
                    action='store_true',
                    help='If in the unbound transformation you are creating the ligand in a box of water use this flag')

# not_generate_csv is a store_false flag
# so it is true when you have to generate the csv files
# and false when you should not
parser.add_argument('--not-generate-csv',
                    action='store_false',
                    help='If you use this flag the work vs lambda CSV files will not be generated, '
                    'generating them tends to be time consuming and they tend to be big files')

parsed_input = parser.parse_args()

bound_dir = Path(parsed_input.bound_dir)
unbound_dir = Path(parsed_input.unbound_dir)

if parsed_input.md_program == 'gromacs':

    def _get_files_for_dir_gromacs(directory, order=['q', 'vdw']):

        output_files = []

        for prefix in order:

            tmp = sorted([i.resolve() for i in directory.glob(f'RESTART/{prefix}*.xvg') if 'pull' not in str(i)])
            
            i = 0
            while (directory / f'Extra_RESTART{i}').exists():
                
                tmp += sorted([i.resolve() for i in directory.glob(f'Extra_RESTART{i}/{prefix}*.xvg') if 'pull' not in str(i)])
                i += 1

            output_files.append(tmp)

            if not tmp:
                return None

        if len(output_files) == 1:
            output_files = output_files[0]
        
        return output_files


    # Bound files
    # There are many possible situations
    possible_filenames = (
        ('q', 'vdw'),
        ('b_1_', 'b_2_'),
        ('b',)
    )

    bound_files = None
    for p in possible_filenames:
        bound_files = _get_files_for_dir_gromacs(bound_dir, order=p)
        if bound_files is not None:
            break
    else:
        raise RuntimeError(f"Could not find the xvg files in {bound_dir}")

    # Unbound files
    # There are many possible situations
    possible_filenames = (
        ('vdw', 'q'),
        ('u_1_', 'u_2_'),
        ('u',)
    )

    unbound_files = None
    for p in possible_filenames:
        unbound_files = _get_files_for_dir_gromacs(unbound_dir, order=p)
        if unbound_files is not None:
            break
    else:
        raise RuntimeError(f"Could not find the xvg files in {unbound_dir}")

    #################################
    # Volume free energy correction
    try:
        volume_correction = gromacs_volume_correction(directory=bound_dir, temperature=parsed_input.temperature)
    
        with open('volume_correction.dat', 'w') as f:
            f.write('# Volume correction in Kcal/mol, must be added to the dissociation '
            'free energy or subtracted from the binding free energy\n'
            f'{volume_correction:.18e}\n')

        del volume_correction

    except RuntimeError:
        with open('volume_correction.dat', 'w') as f:
            f.write('# No file found\n')

    print('created volume_correction.dat with the volume correction')
    ##################################


################################################
# Do vDSSB
if parsed_input.kind_of_process == 'vdssb':

    print('Calculating Jarzinski free energy, will take some time')
    try:
        obj = superclasses.JarzynskiVDSSBPostProcessingPipeline(
            bound_state_dhdl=bound_files,
            unbound_state_dhdl=unbound_files,
            temperature=parsed_input.temperature,
            md_program=parsed_input.md_program,
            bound_creation=parsed_input.bound_creation,
            unbound_creation=parsed_input.unbound_creation)

        free_energy, std = obj.execute()

        print(f'Jarzynski free energy {free_energy:.18e} Kcal/mol\nCI95 {std*1.96:.18e}')
    
    except Exception as e:
        warnings.warn(f"One or more calculations failed because of\n{traceback.format_exc()}")

    print('Calculating Gaussian Mixtures (EM) free energy, will take some time')
    for i in range(3):
        try:
            obj = superclasses.GaussianMixturesVDSSBPostProcessingPipeline(
                bound_state_dhdl=bound_files,
                unbound_state_dhdl=unbound_files,
                temperature=parsed_input.temperature,
                md_program=parsed_input.md_program,
                n_gaussians=i + 1,
                bound_creation=parsed_input.bound_creation,
                unbound_creation=parsed_input.unbound_creation)

            free_energy, std = obj.execute()

            print(f'{i + 1} Gaussian mixtures (EM) free energy {free_energy:.18e} Kcal/mol\nCI95 {std*1.96:.18e}')
        
        except Exception as e:
            warnings.warn(f"One or more calculations failed because of\n{traceback.format_exc()}")

##########################################################
# Do normal FSDAM
else:

    # Jarzynski
    print('Calculating Jarzinski free energy, will take some time')

    try:
        unbound_obj = superclasses.JarzynskiPostProcessingAlchemicalLeg(
            dhdl_files=unbound_files,
            temperature=parsed_input.temperature,
            md_program=parsed_input.md_program,
            creation=parsed_input.unbound_creation)

        unbound_free_energy, unbound_std = unbound_obj.execute()

        shutil.move(f'{str(unbound_obj)}_free_energy.dat',
                    'unbound_' + f'{str(unbound_obj)}_free_energy.dat')

        shutil.move('work_values.dat',
                    'unbound_work_values.dat')

        shutil.move(f'{str(unbound_obj)}_jarzanski_bias.dat', 
            f'unbound_{str(unbound_obj)}_jarzanski_bias.dat')

        print(f'Jarzynski unbound free energy {unbound_free_energy}\n' f'CI95 {1.96*(unbound_std)}')
    

        bound_obj = superclasses.JarzynskiPostProcessingAlchemicalLeg(
            dhdl_files=bound_files,
            temperature=parsed_input.temperature,
            md_program=parsed_input.md_program,
            creation=parsed_input.bound_creation)

        bound_free_energy, bound_std = bound_obj.execute()

        shutil.move(f'{str(bound_obj)}_free_energy.dat',
                    'bound_' + f'{str(bound_obj)}_free_energy.dat')

        shutil.move('work_values.dat',
                    'bound_work_values.dat')

        shutil.move(f'{str(unbound_obj)}_jarzanski_bias.dat', 
            f'bound_{str(unbound_obj)}_jarzanski_bias.dat')

        print(f'Jarzynski bound free energy {bound_free_energy}\n' f'CI95 {1.96*(bound_std)}')

        # Will always print the unbinding free energy
        if parsed_input.bound_creation:
            if parsed_input.unbound_creation:
                total_free_energy = - bound_free_energy + unbound_free_energy
            else:
                total_free_energy = - (bound_free_energy + unbound_free_energy)
        else:
            if parsed_input.unbound_creation:
                total_free_energy = bound_free_energy + unbound_free_energy
            else:
                total_free_energy = bound_free_energy - unbound_free_energy

        total_std = bound_std + unbound_std

        with open(f'{str(bound_obj)}_total_free_energy.dat', 'w') as f:

            f.write(
                '#total unbinding free energy (no volume correction done) in Kcal mol\n'
                '#Dg  STD   CI95%\n'
                f'{total_free_energy:.18e} {total_std:.18e} {1.96*(total_std):.18e}\n'
            )

        print(f'Jarzynski total free energy {total_free_energy}\n' f'CI95 {1.96*(total_std)}')
    
    except Exception as e:
        warnings.warn(f"One or more calculations failed because of\n{traceback.format_exc()}")


    # EM gaussian mixtures
    print('Calculating Gaussian mixtures (EM) free energy, will take some time')
    for i in range(3):
        try:
            unbound_obj = superclasses.GaussianMixturesPostProcessingAlchemicalLeg(
                dhdl_files=unbound_files,
                temperature=parsed_input.temperature,
                md_program=parsed_input.md_program,
                creation=parsed_input.unbound_creation,
                n_gaussians=i + 1)

            unbound_free_energy, unbound_std = unbound_obj.execute()

            shutil.move(f'{str(unbound_obj)}_free_energy.dat',
                        'unbound_' + f'{str(unbound_obj)}_free_energy.dat')

            shutil.move('work_values.dat',
                        'unbound_work_values.dat')

            print(f'{i + 1} Gaussian mixtures (EM) unbound free energy {unbound_free_energy}\n' f'CI95 {1.96*(unbound_std)}')

            bound_obj = superclasses.GaussianMixturesPostProcessingAlchemicalLeg(
                dhdl_files=bound_files,
                temperature=parsed_input.temperature,
                md_program=parsed_input.md_program,
                creation=parsed_input.bound_creation,
                n_gaussians=i + 1)

            bound_free_energy, bound_std = bound_obj.execute()

            shutil.move(f'{str(bound_obj)}_free_energy.dat',
                        'bound_' + f'{str(bound_obj)}_free_energy.dat')

            shutil.move('work_values.dat',
                        'bound_work_values.dat')

            print(f'{i + 1} Gaussian mixtures (EM) bound free energy {bound_free_energy}\n' f'CI95 {1.96*(bound_std)}')

            # Will always print the unbinding free energy
            if parsed_input.bound_creation:
                if parsed_input.unbound_creation:
                    total_free_energy = - bound_free_energy + unbound_free_energy
                else:
                    total_free_energy = - (bound_free_energy + unbound_free_energy)
            else:
                if parsed_input.unbound_creation:
                    total_free_energy = bound_free_energy + unbound_free_energy
                else:
                    total_free_energy = bound_free_energy - unbound_free_energy

            total_std = bound_std + unbound_std

            with open(f'{str(bound_obj)}_total_free_energy.dat', 'w') as f:

                f.write(
                    '#total unbinding free energy (no volume correction done) in Kcal mol\n'
                    '#Dg   STD   CI95%\n'
                    f'{total_free_energy:.18e} {total_std:.18e} {1.96*(total_std):.18e}\n'
                )

            print(f'{i + 1} Gaussian mixtures (EM) total free energy {total_free_energy}\n' f'CI95 {1.96*(total_std)}')

        except Exception as e:
            warnings.warn(f"One or more calculations failed because of\n{traceback.format_exc()}")

##########################################################

# not_generate_csv is a store_false flag
# so it is true when you have to generate the csv files
# and false when you should not
if parsed_input.not_generate_csv:
    print('Creating csv files with work vs lambda (useful for plotting)')
    # Bound
    bound_csv = integrate_works.make_work_vs_lambda_csv(work_files=bound_files,
                                md_program=parsed_input.md_program,
                                creation=parsed_input.bound_creation,
                                num_runs=len(bound_files))

    shutil.move(bound_csv,
                f'bound_{bound_csv}')

    # Unbound
    unbound_csv = integrate_works.make_work_vs_lambda_csv(work_files=unbound_files,
                                md_program=parsed_input.md_program,
                                creation=parsed_input.unbound_creation,
                                num_runs=len(bound_files))

    shutil.move(unbound_csv,
                f'unbound_{bound_csv}')

    print(f'Created {bound_csv} and {unbound_csv}')

# Calculate charge correction

# Get the important info
# and therefore also the ligand resname
unbound_imp_info_file = Path(parsed_input.unbound_dir) / 'important_info.dat'

bound_imp_info_file = Path(parsed_input.bound_dir) / 'important_info.dat'

unbound_imp_info = {}
bound_imp_info = {}

with unbound_imp_info_file.open() as f:
    for line in f:
        line = line.strip()
        if line:
            line = line.split('=')

            unbound_imp_info[line[0].strip()] = line[1].strip()

with bound_imp_info_file.open() as f:
    for line in f:
        line = line.strip()
        if line:
            line = line.split('=')

            bound_imp_info[line[0].strip()] = line[1].strip()

# Using the only solvent box go get the volume, as the volume changes while putting the ligand this is an
# arbitrary approximation 
unbound_vol = mdtraj.load(str(Path(parsed_input.unbound_dir) / unbound_imp_info['only_solvent_gro'])).unitcell_volumes[0]

bound_vol = mdtraj.load(str(Path(parsed_input.bound_dir) / bound_imp_info['gro_file'])).unitcell_volumes[0]

# nm3 to angstrom3
unbound_vol *= 1000
bound_vol *= 1000

unbound_charges = _charge_correction.get_charges_with_parmed(
                                    str(Path(parsed_input.unbound_dir) / unbound_imp_info['top_file']),
                                    xyz=str(Path(parsed_input.unbound_dir) / unbound_imp_info['gro_file']))

unbound_charges = sum(unbound_charges)

bound_charges = _charge_correction.get_charges_with_parmed(
                                    str(Path(parsed_input.bound_dir) / bound_imp_info['top_file']),
                                    xyz=str(Path(parsed_input.bound_dir) / bound_imp_info['gro_file']))

bound_charges = sum(bound_charges)

homogeneus_correction = _charge_correction.homogeneus_charge_correction_vDSSB(
        host_charge=bound_charges - unbound_charges,
        guest_charge=unbound_charges,
        host_guest_box_volume=bound_vol,
        only_guest_box_volume=unbound_vol)

# hartee to kcal/mol
homogeneus_correction *= 627.5

globular_correction = _charge_correction.globular_protein_correction(
        pdb_file=str(Path(parsed_input.bound_dir) / bound_imp_info['gro_file']),
        host_charge=bound_charges - unbound_charges,
        guest_charge=unbound_charges,
        ligand=f'resname {unbound_imp_info["ligand_resname"]}')

lines = ('# This charge corrections are in Kcal/mol and is for the vDSSB method\n'
    '# They should be added to the dissociacion free energy or removed from the binding free energy\n'
    '# The calculation took for granted that there are no counterions in the unbound system\n'
    '# The homogeneus correction should always be used in case of charged ligands\n'
    '# The globular correction is only valid with globular proteins (it assumes the protein as a perfact sphere)\n'
    '# All the charges are in atomic charges\n'
    '# Check ref https://doi.org/10.1063/5.0086640 and https://pubs.acs.org/doi/abs/10.1021/ct400626b \n'
    '\n'
    'ATTENTION! In gromacs the homogeneus correction is done under the hood by gromacs itself\n'
    'See the function ewald_charge_correction in the source file src/gromacs/ewald/ewald.cpp of the gromacs MD program\n'
    '\n'
    f'ligand_charge = {unbound_charges:.18e}\n'
    f'host_charge = {bound_charges - unbound_charges:.18e}\n'
    f'homogeneus_correction = {homogeneus_correction:.18e}\n'
    f'globular_correction = {globular_correction:.18e}\n')

with open('charge_correction.dat', 'w') as f:
    f.write(lines)

print('Created charge_correction.dat')
