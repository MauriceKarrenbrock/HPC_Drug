######################################################################################
# Copyright (c) 2021-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

import argparse
from openmmforcefields.generators import SMIRNOFFTemplateGenerator, GAFFTemplateGenerator

from HPC_Drug.MD.parametrize_with_openmm import parametrize_protein_and_ligand

parser = argparse.ArgumentParser(
    description='This script parametrizes a given system with openmmforcefields/openmm '
    'to work it needs a pdb with the (solvated) protein, a sdf file of the protonated ligand '
    'and a box of water (optional)',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--protein',
    action = "store",
    required=True,
    type=str,
    help = "A pdb file containing the protein, water and ions but not the ligand")


parser.add_argument('--ligand',
    action = "store",
    required=True,
    type=str,
    help = "A sdf file of the protonated ligand, if it doesn't have the right coordinates to "
    "bind the protein the parametrization will still be succesful but the obtained coordinate files will "
    "of course be useless.\n"
    "If you give the ligand in another format it will be converted to sdf with openbabel at your own risk")

parser.add_argument('--ligand-resname',
    action = "store",
    default="LIG",
    type=str,
    help = "The residue name of the ligand (key sensitive)")

parser.add_argument('--ligand-protonation-ph',
    action = "store",
    default=None,
    type=float,
    help = "If the ligand sdf file is not protonated you can give the pH at which "
    "the program should protonate it with openbabel, per default the ligand remains "
    "untouched")

parser.add_argument('--forcefield-files',
    action = "store",
    default='amber/ff14SB.xml,amber/tip3p_standard.xml',
    type=str,
    help = "The forcefield files in a comma separated list for the protein, solvent and ions in xml format "
    "compatible with simtk.openmm.Forcefield they can be the available ones in openmm and openmmforcefields "
    "or a custom one in the working directory. There is no limit to the number of files. "
    "Do not give the ligand forcefield here")

parser.add_argument('--ligand-forcefield',
    action = "store",
    default=None,
    type=str,
    help = "The forcefield for the ligand, can be any of "
    "the supported ones by openmmforcefields.generators, at the moment Gaff and the "
    "openff initiative ones (es SMIRNOFF) are the only ones supported\n"
    "The default is the newest Gaff forcefield avaiable\n"
    "The supported ones with your openmmforcefields version are \n\n"
    f"{' '.join(GAFFTemplateGenerator.INSTALLED_FORCEFIELDS)}\n\n"
    f"{' '.join(SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS)}\n\n")

parser.add_argument('--water-box-output',
    action = "store",
    default='only_water',
    type=str,
    help = "The prefix name for the files of a box of only water that will be generated "
    "pre existing files will be overwritten!")

parser.add_argument('--water-type',
    action = "store",
    default='tip3p',
    type=str,
    help = "The water model to use when creating the water box, only water types supported "
    "by openmmtools.testsystems.WaterBox")

parser.add_argument('--water-box-edge',
    action = "store",
    default=3.0,
    type=float,
    help = "The box edge of the water box in nm, default 3.0 (optional)")

parser.add_argument('--complex-output',
    action = "store",
    default='complex',
    type=str,
    help = "the prefix for all the outputs concerning the complex "
    "es complex.pdb, complex.top, ..."
    "attention not to overwrite existing files, no check will be done!")

parser.add_argument('--ligand-output',
    action = "store",
    default='LIG',
    type=str,
    help = "the prefix for all the outputs concerning the ligand "
    "es LIG.pdb, LIG.top, ..."
    "attention not to overwrite existing files, no check will be done!")


parsed_input = parser.parse_args()

parsed_input.forcefield_files = parsed_input.forcefield_files.split(',')

parametrize_protein_and_ligand(protein_water=parsed_input.protein,
        ligand=parsed_input.ligand,
        forcefields=parsed_input.forcefield_files,
        ligand_forcefield=parsed_input.ligand_forcefield,
        complex_prefix=parsed_input.complex_output,
        ligand_prefix=parsed_input.ligand_output,
        ligand_resname=parsed_input.ligand_resname,
        water_prefix=parsed_input.water_box_output,
        box_edge=parsed_input.water_box_edge,
        water_model=parsed_input.water_type,
        ligand_ph=parsed_input.ligand_protonation_ph)
