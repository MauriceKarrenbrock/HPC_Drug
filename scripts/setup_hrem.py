######################################################################################
# Copyright (c) 2021-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

import argparse
import parmed
from HPC_Drug.MD.gromacs.hrem import ComplexHREMInput, LigandHREMInput

parser = argparse.ArgumentParser(
    description='This script will set up an hamiltonian replica exchange both for '
    'the protein ligand system and for the ligand alone. '
    'THE INPUT FILES GENERATED ARE VERY GENERIC AND SHALL NEVER BE USED BLINDLY!',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--temperature',
    action = 'store',
    default=298.15,
    type=float,
    help = 'Temperature in Kelvin')

parser.add_argument('--md-program',
    action = 'store',
    default='gromacs',
    choices = ['gromacs'],
    type=str,
    help = 'For which MD prgram the HREM has to be set up. '
    'In case of Groamcs a heavy dummy atom will be put on the center of mass of the protein '
    'because of some PBC problems')

parser.add_argument('--complex-coordinates',
    action = 'store',
    default='',
    type=str,
    help = 'A pdb or gro file '
    'that contains the solvated protein ligand complex')

parser.add_argument('--ligand-coordinates',
    action = 'store',
    default='',
    type=str,
    help = 'A pdb or gro file '
    'that contains the ligand in vacuum in a big box (es 15nm)')

parser.add_argument('--water-coordinates',
    action = 'store',
    default='',
    type=str,
    help = 'A pdb or gro file '
    'that contains a box of water')

parser.add_argument('--complex-topology',
    action = 'store',
    default='',
    type=str,
    help = 'A a gromacs top file '
    'that contains the solvated protein ligand complex')

parser.add_argument('--ligand-topology',
    action = 'store',
    default='',
    type=str,
    help = 'A a gromacs top file '
    'that contains the ligand in vacuum')

parser.add_argument('--water-ligand-topology',
    action = 'store',
    default='',
    type=str,
    help = 'A a gromacs top file '
    'for the box of water and the ligand (in this order)')

parser.add_argument('--ligand-resname',
    action = 'store',
    required=True,
    type=str,
    help = 'The residue name of the ligand (case sensitive)')

parser.add_argument('--ligand-resnumber',
    action = 'store',
    default=-100,
    type=int,
    help = 'The residue number of the ligand as written in the complex pdb (resSeq)')

parser.add_argument('--complex-batteries',
    action = 'store',
    default=-1,
    type=int,
    help = 'How many separate HREM you want to run in parallel, '
    'if it is -1 it will be chosen automatically')

parser.add_argument('--ligand-batteries',
    action = 'store',
    default=-1,
    type=int,
    help = 'How many separate HREM you want to run in parallel, '
    'if it is -1 it will be chosen automatically')

parser.add_argument('--complex-timestep',
    action = 'store',
    default=0.002,
    type=float,
    help = 'Timestep for the HREM')

parser.add_argument('--ligand-timestep',
    action = 'store',
    default=0.002,
    type=float,
    help = 'Timestep for the HREM')

parser.add_argument('--complex-nsteps',
    action = 'store',
    default=-1,
    type=int,
    help = 'Number of steps for the HREM divided among the number of batteries '
    'if it is -1 it will be chosen automatically to get circa 32 ns in total')

parser.add_argument('--ligand-nsteps',
    action = 'store',
    default=-1,
    type=int,
    help = 'Number of steps for the HREM divided among the number of batteries '
    'if it is -1 it will be chosen automatically to get circa 32 ns in total')

parser.add_argument('--complex-nreplicas',
    action = 'store',
    default=8,
    type=int,
    help = 'Number of replicas for the HREM')

parser.add_argument('--ligand-nreplicas',
    action = 'store',
    default=8,
    type=int,
    help = 'Number of replicas for the HREM')

parser.add_argument('--no-preprocess-topologies',
    action = 'store_false',
    help = 'If you put this flag your topologies will not be preprocessed by gromacs in order '
    'to remove #import statements (it means that there are none already)')

parser.add_argument('--complex-constraints',
    action = 'store',
    default='h-bonds',
    choices = ['h-bonds', 'all-bonds', 'none'],
    type=str,
    help = 'The constraints to use for the complex HREM '
    'the sintax is the same as in gromacs')

parser.add_argument('--ligand-constraints',
    action = 'store',
    default='h-bonds',
    choices = ['h-bonds', 'all-bonds', 'none'],
    type=str,
    help = 'The constraints to use for the ligand HREM '
    'the sintax is the same as in gromacs')

parser.add_argument('--complex-scaling-top-basis',
    action = 'store',
    default=0.2,
    type=float,
    help = 'The basis for the geometrical progression used to scale the topologies '
    'base**(i/n_replicas-1) for i in range(i)')

parser.add_argument('--ligand-scaling-top-basis',
    action = 'store',
    default=0.1,
    type=float,
    help = 'The basis for the geometrical progression used to scale the topologies '
    'base**(i/n_replicas-1) for i in range(i)')

parser.add_argument('--cores-per-node',
    action = 'store',
    default=64,
    type=int,
    help = 'Number of CPUs per node')

parser.add_argument('--gpus-per-node',
    action = 'store',
    default=1,
    type=int,
    help = 'Number of GPUs per node')

parsed_input = parser.parse_args()

if parsed_input.complex_batteries == -1:
    parsed_input.complex_batteries = None

if parsed_input.ligand_batteries == -1:
    parsed_input.ligand_batteries = None

if parsed_input.complex_nsteps == -1:
    parsed_input.complex_nsteps = None

if parsed_input.ligand_nsteps == -1:
    parsed_input.ligand_nsteps = None

if parsed_input.md_program == 'gromacs':

    if parsed_input.complex_topology and parsed_input.complex_coordinates and (parsed_input.ligand_resnumber != -100):
        # Convert coordinates to gro
        # Complex
        pmd = parmed.load_file(parsed_input.complex_topology, xyz=parsed_input.complex_coordinates)
        pmd.save(parsed_input.complex_coordinates + '.gro', overwrite=True)
        parsed_input.complex_coordinates = parsed_input.complex_coordinates + '.gro'

        _, _, complex_dir = ComplexHREMInput(
            pdb_file=parsed_input.complex_coordinates,
            top_file=parsed_input.complex_topology,
            ligand_resname=parsed_input.ligand_resname,
            ligand_resnumber=parsed_input.ligand_resnumber,
            number_of_cores_per_node=parsed_input.cores_per_node,
            gpus_per_node=parsed_input.gpus_per_node,
            number_of_replicas=parsed_input.complex_nreplicas,
            batteries=parsed_input.complex_batteries,
            n_steps=parsed_input.complex_nsteps,
            timestep=parsed_input.complex_timestep,
            constraints=parsed_input.complex_constraints,
            temperature=parsed_input.temperature,
            top_scaling_basis=parsed_input.complex_scaling_top_basis,
            preprocess_topology=parsed_input.no_preprocess_topologies
        ).execute()
    
    else:
        print("Complex HREM will not be set up because you did not input all the needed parameters")
        complex_dir = ''


    if parsed_input.ligand_topology and parsed_input.ligand_coordinates and parsed_input.water_coordinates and parsed_input.water_ligand_topology:
        # Convert coordinates to gro
        # Ligand
        pmd = parmed.load_file(parsed_input.ligand_topology, xyz=parsed_input.ligand_coordinates)
        pmd.save(parsed_input.ligand_coordinates + '.gro', overwrite=True)
        parsed_input.ligand_coordinates = parsed_input.ligand_coordinates + '.gro'
        
        ligand_dir = LigandHREMInput(
            ligand_gro_file=parsed_input.ligand_coordinates,
            ligand_top_file=parsed_input.ligand_topology,
            water_gro_file=parsed_input.water_coordinates,
            water_ligand_top_file=parsed_input.water_ligand_topology,
            ligand_resname=parsed_input.ligand_resname,
            number_of_cores_per_node=parsed_input.cores_per_node,
            gpus_per_node=parsed_input.gpus_per_node,
            number_of_replicas=parsed_input.ligand_nreplicas,
            batteries=parsed_input.ligand_batteries,
            n_steps=parsed_input.ligand_nsteps,
            timestep=parsed_input.ligand_timestep,
            constraints=parsed_input.ligand_constraints,
            temperature=parsed_input.temperature,
            top_scaling_basis=parsed_input.ligand_scaling_top_basis,
            preprocess_topology=parsed_input.no_preprocess_topologies
        ).execute()
    
    else:
        print("Ligand HREM will not be set up because you did not input all the needed parameters")
        ligand_dir = ''

    print(f'Two directories were created\n{complex_dir}\n{ligand_dir}')