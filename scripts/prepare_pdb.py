######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

import argparse
import shutil
import warnings

from simtk import unit
from simtk.openmm.app import Modeller, PDBFile, ForceField
import mdtraj
import numpy as np
from openff.toolkit.topology import Molecule
import pdbfixer

from PythonPDBStructures import geometry as _geo

from HPC_Drug.PDB.repair_pdb.repair import repair
from HPC_Drug.PDB.select_model_chain import select_model_chain
from HPC_Drug.PDB.add_chain_id import add_chain_id
from HPC_Drug.PDB.remove_disordered_atoms import remove_disordered_atoms
from HPC_Drug.PDB.structural_information.scan_structure import get_organic_ligands_with_no_header
from HPC_Drug.MD import openbabel_utils
from HPC_Drug.PDB.structural_information.scan_structure import get_metal_binding_residues_with_no_header, get_disulf_bonds_with_no_header
from HPC_Drug.PDB.complex_selections import remove_trash
from HPC_Drug.MD import residue_renaming
from HPC_Drug.PDB import biopython as _bio_utils

parser = argparse.ArgumentParser(
    description='This script helps you clean up and prepare a protein-ligand '
    'pdb file for molecular dynamics, it can do many things for you like '
    'repairing missing atoms/residues with pdbfixer, select a certain model and chain '
    'from a bigger pdb file, add a chain id if your pdb is missing it (it is needed for many parsers), '
    'rotate the system in the reference frame of the moment of inertia tensor, '
    'remove disordered atoms (always done),  etc. '
    'but pay attention that the residue renaming and subsequent protonation are done with '
    'quick and dirty defaults so if you know you have a more complex system use other tools',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

parser.add_argument('--input-pdb',
    action = 'store',
    required=True,
    type=str,
    help = 'The input PDB file, attention: all the hetero atoms must be labelled HETATM '
    'and all the protein atoms ATOM')

parser.add_argument('--output-pdb',
    action = 'store',
    required=True,
    type=str,
    help = 'The output PDB file, must be different than the input one')

parser.add_argument('--select-chain',
    action = 'store',
    default=None,
    type=str,
    help = 'The chain id to select, case sensitive, default=no selection done')

parser.add_argument('--select-model',
    action = 'store',
    default=0,
    type=int,
    help = 'The model id to select, default=The first one (0)')

parser.add_argument('--add-chain-id',
    action = 'store',
    default=None,
    type=str,
    help = 'The chain id to add to the pdb file, case sensitive, default=no chain id added')

parser.add_argument('--repair-pdb',
    action = 'store_true',
    help = 'If to add missing atoms/residues, and change non standard residue names '
    'in standard ones with pdbfixer, it is a flag')

parser.add_argument('--remove-trash',
    action = 'store_true',
    help = 'If to check for ions that are usually there for '
    'the crystallization process and remove them (of course check the output), it is a flag')

parser.add_argument('--create-ligand-sdf',
    action = 'store_true',
    help = 'DEPRECATED, if you do not have the SDF files for the ligand(s) in the PDB '
    'this program can generate them with openbabel from the PDB file, at your own risk. '
    'If you got the PDB from the databank you can download the sdf files with the '
    'right coordinates from there')

parser.add_argument('--rotate-system',
    action = 'store_true',
    help = 'If to rotate the protein-ligand system in its moment of inertia tensor '
    'reference system, useful for solvating in a smaller box')

parser.add_argument('--sdf-files-to-rotate',
    action = 'store',
    default=None,
    type=str,
    help = 'The sdf files to rotate (comma separated list) if the --rotate-system flag is used, '
    'in this way their coordinates will still be valid, if you have put the '
    '--create-ligand-sdf flag thos sdf files will be rotated automatically')

parser.add_argument('--residue-renaming',
    action = 'store',
    default=None,
    choices = ['amber'],
    type=str,
    help = 'According to your forcefield choice residues binding metallic ions are renamed, '
    'and some residues will also be renamed according to the pH, default=no renaming')

parser.add_argument('--ph',
    action = 'store',
    default=7.0,
    type=float,
    help = 'The pH value')

parser.add_argument('--add-hydrogens',
    action = 'store_true',
    help = 'If to add hydrogens with openmm.app.Modeller.addHydrogens, '
    'default=do not add hydrogens')

parser.add_argument('--custom-hydrogen-definitions',
    action = 'store',
    default=None,
    type=str,
    help = 'Files to give to openmm.app.Modeller with the loadHydrogenDefinitions method, '
    'comma separated list')

parser.add_argument('--solvate',
    action = 'store_true',
    help = 'If to solvate the system in a cubic box of water')

parser.add_argument('--solvate-padding',
    action = 'store',
    default=1.2,
    type=float,
    help = 'Padding in nm to apply for the solvation')

parser.add_argument('--water-type',
    action = 'store',
    default='tip3p',
    type=str,
    help = 'The kind of water to add, it must be supported by openmm, case sensitive')

parser.add_argument('--neutralize',
    action = 'store_true',
    help = 'If you want the water box to be neutralized, it is a flag')

parser.add_argument('--forcefield-files',
    action = 'store',
    default='amber/ff14SB.xml,amber/tip3p_standard.xml',
    type=str,
    help = 'The forcefield files in a comma separated list for the protein, solvent and ions in xml format '
    'compatible with simtk.openmm.Forcefield they can be the available ones in openmm and openmmforcefields '
    'or a custom one in the working directory. There is no limit to the number of files. '
    'Do not give the ligand forcefield here. '
    'The forcefield is only used to add water and hydrogens, so if you use another one later '
    'it should not be a big deal if you do an energy minimization')

parsed_input = parser.parse_args()

if parsed_input.sdf_files_to_rotate:
    parsed_input.sdf_files_to_rotate = parsed_input.sdf_files_to_rotate.split(',')
else:
    parsed_input.sdf_files_to_rotate = []

if parsed_input.custom_hydrogen_definitions:
    parsed_input.custom_hydrogen_definitions = parsed_input.custom_hydrogen_definitions.split(',')

parsed_input.forcefield_files = parsed_input.forcefield_files.split(',')

# Copy input PDB in output PDB in order to be sure not to modify it later on
shutil.copy(parsed_input.input_pdb, parsed_input.output_pdb)

output_pdb = parsed_input.output_pdb

# Download
# TODO, maybe

# Repair
if parsed_input.repair_pdb:
    repair(input_file=output_pdb,
        output_file=output_pdb,
        repairing_method="pdbfixer")

# select chain and model
select_model_chain(input_pdb=output_pdb,
                    output_pdb=output_pdb,
                    model=parsed_input.select_model,
                    chain=parsed_input.select_chain)

# add chain id
if parsed_input.add_chain_id:
    add_chain_id(pdb_file=output_pdb, chain=parsed_input.add_chain_id)


# remove disordered atoms
remove_disordered_atoms(input_pdb=output_pdb,
                    output_pdb=output_pdb)
                    
###########################################################
# PATCH
###########################################################

#a patch to add possible missing TER lines in pdb files
#adressing issues #2733 and #2736 on biopython github
#Will be removed when the issues will be solved
with open(output_pdb, 'r') as f:
    lines = f.readlines()

with open(output_pdb, 'w') as f:
    for i in range(len(lines)-1):
        if lines[i][0:4] == 'ATOM' and lines[i + 1][0:6] == 'HETATM':
            #writes a TER line on the last ATOM of each chain (not robust for not natural AA)
            lines[i] = lines[i].strip() + '\n' + "{0:<6}{1}{2}{3}".format("TER", lines[i][6:11], ' '*6, lines[i][17:26])
        
        f.write(f"{lines[i].strip()}\n")

    f.write(f"{lines[len(lines)-1].strip()}\n")

###########################################################
# END PATCH
###########################################################

# extract ligands
extracted_ligands = []
if parsed_input.create_ligand_sdf:
    lig_names_resSeq = get_organic_ligands_with_no_header(pdb_file=output_pdb, trash=[])

    for lig in lig_names_resSeq:
        lig_name_string = f'{lig[0]}_resSeq{lig[1]}'
        traj =  mdtraj.load(output_pdb)
        traj = traj.atom_slice(traj.top.select(f'resSeq {lig[1]}'))
        traj.save(lig_name_string + '.pdb', force_overwrite=True)

        openbabel_utils.convert_and_protonate_file_to_sdf(input_file=lig_name_string + '.pdb',
            sdf_file=lig_name_string + '.sdf',
            ph=None,
            ligand_resname=lig[0])

        extracted_ligands.append(lig_name_string + '.sdf')


# remove trash
if parsed_input.remove_trash:
    remove_trash(input_file=output_pdb,
        output_file=output_pdb)


# rotate prot and lig-sdf
if parsed_input.rotate_system:
    rot_matrix = _geo.get_inertia_eigenvectors(output_pdb)
    rot_matrix = np.linalg.inv(rot_matrix)

    traj = mdtraj.load(output_pdb)

    # rotate pdb
    traj.xyz[0] = _geo.rotate_coordinates(coordinates=traj.xyz[0],
                                        rot_matrix=rot_matrix,
                                        check_reflections=True)

    traj.save(output_pdb, force_overwrite=True)

    for sdf in (extracted_ligands + parsed_input.sdf_files_to_rotate):

        # Sometimes openff-toolkit doen't like openbabel's generated sdf files
        try:
            lig = Molecule.from_file(sdf, 'sdf', allow_undefined_stereo=True)
            # rotate sdf (if given)
            lig.conformers[0] = _geo.rotate_coordinates(coordinates=lig.conformers[0] / unit.nanometers,
                                                rot_matrix=rot_matrix,
                                                check_reflections=True) * unit.nanometers

            # Attention, sometimes they result translated for some reason
            lig.to_file('rotated_' + sdf, 'sdf')

        except Exception as e:
            warnings.warn(f'Could not rotate SDF file {sdf} because of: ' + str(e))


# Rename
if parsed_input.residue_renaming:
    structure = _bio_utils.get_structure(output_pdb)

    substitutions_dict = get_metal_binding_residues_with_no_header(structure,
                                                protein_chain = None,
                                                protein_model = None)

    substitutions_dict, sulf_bonds = get_disulf_bonds_with_no_header(structure,
                                                protein_chain = None,
                                                protein_model = None)


    structure = residue_renaming.ResidueRenamer(structure, #biopython structure
                    substitutions_dict,
                    forcefield=parsed_input.residue_renaming,
                    substitution = "standard",
                    ph = parsed_input.ph).execute()
    
    structure = None

# protonate
structure = PDBFile(output_pdb)
modeller = Modeller(structure.topology, structure.positions)

if parsed_input.add_hydrogens:

    if parsed_input.custom_hydrogen_definitions:
        modeller.loadHydrogenDefinitions(parsed_input.custom_hydrogen_definitions)

    variants = []

    variants_dict = dict(CYM='CYX')

    for residue in structure.topology.residues():
        variants.append(variants_dict.get(residue.name, None))

    modeller.addHydrogens(pH=parsed_input.ph, variants=variants)

    # Write it down in case solvation fails
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, file=f)

# solvate
if parsed_input.solvate:

    fixer = pdbfixer.PDBFixer(filename=output_pdb)

    # This is a small hack in order to have a fast and easy dummy FF
    # In order to add waters
    dummy_forcefield = fixer._createForceField(fixer.topology, True)

    modeller.addSolvent(forcefield=dummy_forcefield,
                        model=parsed_input.water_type,
                        padding=parsed_input.solvate_padding * unit.nanometers,
                        neutralize=parsed_input.neutralize)

    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, file=f)

print('Setup done!')
