######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""Functions to parametrize a system with openmm/openmmforcefields

`parametrize_protein_and_ligand` is a high level function that will take care
about the protein and the ligand
"""

import numpy as np

from simtk import unit
from simtk.openmm.app import PDBFile, ForceField, PME, Modeller
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator, SMIRNOFFTemplateGenerator
from openmmtools import testsystems

import parmed

import HPC_Drug.MD.openbabel_utils as _obabel

def get_ligand_ff_generator(ligand_forcefield, ligand):
    """Returns the right ligand forcefield generator
    
    Parameters
    -------------
    ligand_forcefield : str
        the forcefield to use, if None the latest gaff will be used
    ligand : openff.toolkit.topology.Molecule

    Returns
    -----------
    openmmforcefields.generators.GAFFTemplateGenerator or openmmforcefields.generators.SMIRNOFFTemplateGenerator
    """
    if ligand_forcefield is  None:
        ligand_ff_generator = GAFFTemplateGenerator(molecules=ligand)

    elif ligand_forcefield in GAFFTemplateGenerator.INSTALLED_FORCEFIELDS:
        ligand_ff_generator = GAFFTemplateGenerator(molecules=ligand, forcefield=ligand_forcefield)

    elif ligand_forcefield in SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS:
        ligand_ff_generator = SMIRNOFFTemplateGenerator(molecules=ligand, forcefield=ligand_forcefield)

    else:
        raise ValueError(f'{ligand_forcefield} is not a supported ligand forcefield\n'
        'the supported ones are:\n'
        f"{' '.join(GAFFTemplateGenerator.INSTALLED_FORCEFIELDS)}\n"
        f"{' '.join(SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS)}")

    return ligand_ff_generator

def create_water_box(box_edge=3.0*unit.nanometers,
    water_model='tip3p'):
    """Creates a water box with openmmtools.testsystems.WaterBox

    Parameters
    -------------
    box_edge : openmm.Quantity or float, default=3.0*unit.nanometers
        if a float is given it will be considered in nm
    water_model : str, default=tip3p
        a water model supported by openmmtools.testsystems.WaterBox
    
    Returns
    -----------
    a box of only water with attributes `positions` and `topology` in openmm format
    """
    if not unit.is_quantity(box_edge):
        box_edge = box_edge * unit.nanometers

    water_box = testsystems.WaterBox(box_edge=box_edge, model=water_model,
                                 constrained=False, nonbondedMethod=PME, ionic_strength=0*unit.molar)

    return water_box

def parametrize_complex(protein_water,
    ligand,
    forcefields=None,
    ligand_forcefield=None,
    output_prefix='complex'):
    """"Parametrizes the protein ligand water ions complex
    and creates .top .pdb .gro .prmtop .inpcrd files
    
    Paramenters
    -------------
    protein_water : str
        the pdb file of the protein, water and ions (not the ligand)
    ligand : str or openff.toolkit.topology.Molecule
        the SDF file of the ligand or a Molecule instance
        if the coordinates of the ligand are not the right ones
        the parametrization will still be succesful but of course the
        obtained structure files will be useless
    forcefields : iterable(str), default=['amber/ff14SB.xml', 'amber/tip3p_standard.xml']
        The forcefield files for the protein, solvent and ions in xml format
        compatible with simtk.openmm.Forcefield they can be the available ones
        in openmm and openmmforcefields or a custom one in the working directory.
        There is no limit to the number of files.
        Do not give the ligand forcefield here
    ligand_forcefield : str, default=the newest gaff FF available
        The forcefield for the ligand, can be any of
        the supported ones by openmmforcefields.generators, at the moment Gaff and the
        openff initiative ones (es SMIRNOFF) are the only ones supported
    output_prefix : str, default='complex'
        the prefix for all the created files, WILL OVERWRITE PRE EXISTSING ONES!
    """

    if not isinstance(ligand, Molecule):
        ligand = Molecule.from_file(ligand)

    if forcefields is None:
        forcefields = ['amber/ff14SB.xml', 'amber/tip3p_standard.xml']

    ligand_ff_generator = get_ligand_ff_generator(ligand_forcefield=ligand_forcefield,
        ligand=ligand)

    system_kwargs = {'constraints': None, 'rigidWater': False, 'nonbondedMethod': PME}

    forcefield = ForceField(*forcefields)

    forcefield.registerTemplateGenerator(ligand_ff_generator.generator)

    protein_water_pdb = PDBFile(protein_water)
    protein_ligand_modeller = Modeller(protein_water_pdb.topology,
        protein_water_pdb.positions)

    protein_ligand_modeller.add(ligand.to_topology().to_openmm(), ligand.conformers[0])

    # Make the protein ligand .top file
    protein_ligand_system = forcefield.createSystem(
        protein_ligand_modeller.topology, **system_kwargs)

    pmd_structure = parmed.openmm.load_topology(protein_ligand_modeller.topology,
    system=protein_ligand_system, xyz=protein_ligand_modeller.positions)

    # Gromacs files
    pmd_structure.save(output_prefix + '.top', overwrite=True)
    pmd_structure.save(output_prefix + '.pdb', overwrite=True)
    pmd_structure.save(output_prefix + '.gro', overwrite=True)

    # Amber files
    pmd_structure.save(output_prefix + '.prmtop', overwrite=True)
    pmd_structure.save(output_prefix + '.inpcrd', overwrite=True)

    # The parmed top files have system1 instead of Protein and when
    # defining a restraint gromacs will complain
    with open(output_prefix + '.top', 'r') as f:
        protein_top_str = f.read()

    with open(output_prefix + '.top', 'w') as f:
        f.write(protein_top_str.replace('system1', 'Protein'))


def parametrize_ligand(ligand,
    water_box,
    water_forcefields=None,
    ligand_forcefield=None,
    ligand_prefix='LIG',
    water_prefix='only_water'):
    """Parametrizes the ligand and a box of water

    Creates gromacs and amber files for the ligand, the water, and water + ligand
    (in this order)
    .top .pdb .gro .prmtop .inpcrd

    Paramenters
    -------------
    ligand : str or openff.toolkit.topology.Molecule
        the SDF file of the ligand or a Molecule instance
    water_box : object with .positions and .topology in openmm format
    water_forcefields : iterable(str), default=['amber/tip3p_standard.xml']
        The forcefield files for the protein, solvent and ions in xml format
        compatible with simtk.openmm.Forcefield they can be the available ones
        in openmm and openmmforcefields or a custom one in the working directory.
        There is no limit to the number of files.
        Do not give the ligand forcefield here
    ligand_forcefield : str, default=the newest gaff FF available
        The forcefield for the ligand, can be any of
        the supported ones by openmmforcefields.generators, at the moment Gaff and the
        openff initiative ones (es SMIRNOFF) are the only ones supported
    ligand_prefix : str, default='LIG'
        the prefix for all the created ligand files, WILL OVERWRITE PRE EXISTSING ONES!
    water_prefix : str, default='LIG'
        the prefix for all the created water files, WILL OVERWRITE PRE EXISTSING ONES!
    """
    if not isinstance(ligand, Molecule):
        ligand = Molecule.from_file(ligand)

    if water_forcefields is None:
        water_forcefields = ['amber/tip3p_standard.xml']
    
    ligand_ff_generator = get_ligand_ff_generator(ligand_forcefield=ligand_forcefield,
        ligand=ligand)

    system_kwargs = {'constraints': None, 'rigidWater': False, 'nonbondedMethod': PME}

    forcefield = ForceField(*water_forcefields)

    forcefield.registerTemplateGenerator(ligand_ff_generator.generator)

    # Only water
    water_box.system = forcefield.createSystem(
        water_box.topology, **system_kwargs)

    pmd_structure = parmed.openmm.load_topology(water_box.topology,
    system=water_box.system, xyz=water_box.positions)

    pmd_structure.save(water_prefix + '.top', overwrite=True)
    pmd_structure.save(water_prefix + '.pdb', overwrite=True)
    pmd_structure.save(water_prefix + '.gro', overwrite=True)

    pmd_structure.save(water_prefix + '.prmtop', overwrite=True)
    pmd_structure.save(water_prefix + '.inpcrd', overwrite=True)

    # Only ligand
    ligand_positions, ligand_topology =  ligand.conformers[0], ligand.to_topology().to_openmm()

    ##########
    # Bring the minimum coordinates to 0 0 0
    tmp_coordinates = np.array(ligand_positions / unit.nanometers)

    for i in range(3):

        tmp_coordinates[:, i] = tmp_coordinates[:, i] - np.amin(tmp_coordinates[:, i])

    ligand_positions = tmp_coordinates *  unit.nanometers

    del tmp_coordinates
    ############

    ligand_topology.setUnitCellDimensions([15.0, 15.0, 15.0] * unit.nanometers)

    ligand_system = forcefield.createSystem(ligand_topology, **system_kwargs)

    pmd_structure = parmed.openmm.load_topology(ligand_topology,
    system=ligand_system, xyz=ligand_positions)

    pmd_structure.save(ligand_prefix + '.top', overwrite=True)
    pmd_structure.save(ligand_prefix + '.pdb', overwrite=True)
    pmd_structure.save(ligand_prefix + '.gro', overwrite=True)

    pmd_structure.save(ligand_prefix + '.prmtop', overwrite=True)
    pmd_structure.save(ligand_prefix + '.inpcrd', overwrite=True)

    # Water ligand topology (water forst ligand second)
    ligand_modeller = Modeller(water_box.topology, water_box.positions)
    ligand_modeller.add(ligand_topology, ligand_positions)

    ligand_system = forcefield.createSystem(
        ligand_modeller.topology, **system_kwargs)

    pmd_structure = parmed.openmm.load_topology(ligand_modeller.topology,
        system=ligand_system, xyz=ligand_modeller.positions)

    pmd_structure.save(water_prefix + '_and_' + ligand_prefix + '.top', overwrite=True)
    pmd_structure.save(water_prefix + '_and_' + ligand_prefix + '.prmtop', overwrite=True)



def parametrize_protein_and_ligand(protein_water,
        ligand,
        forcefields=None,
        ligand_forcefield=None,
        complex_prefix='complex',
        ligand_prefix='LIG',
        ligand_resname='LIG',
        water_prefix='only_water',
        box_edge=3.0*unit.nanometers,
        water_model='tip3p',
        ligand_ph=None):

    _obabel.convert_and_protonate_file_to_sdf(ligand,
        ligand_prefix + '.sdf', ph=ligand_ph, ligand_resname=ligand_resname)

    ligand = Molecule.from_file(ligand_prefix + '.sdf')
    
    water_box = create_water_box(box_edge=box_edge,
        water_model=water_model)

    parametrize_ligand(ligand=ligand,
        water_box=water_box,
        water_forcefields=forcefields,
        ligand_forcefield=ligand_forcefield,
        ligand_prefix=ligand_prefix,
        water_prefix=water_prefix)


    parametrize_complex(protein_water=protein_water,
        ligand=ligand,
        forcefields=forcefields,
        ligand_forcefield=ligand_forcefield,
        output_prefix=complex_prefix)
