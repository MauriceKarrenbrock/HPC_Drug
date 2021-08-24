######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""Functions to parametrize a system with openmm/openmmforcefields
"""

import numpy as np

import pdbfixer
from simtk import unit
from simtk.openmm.app import PDBFile, ForceField, PME, Modeller
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator, SMIRNOFFTemplateGenerator
from openmmtools import testsystems

import parmed
import mdtraj

import HPC_Drug.PDB.complex_selections as _selections
import HPC_Drugs.PDB.openbabel_utils as _obabel

def _scale_positions_and_update_box_vectors(positions, topology):

    coordinates = np.array(positions / unit.nanometers)

    max_coord = []
    for i in range(3):

        coordinates[:, i] = coordinates[:, i] - np.amin(coordinates[:, i])

        max_coord.append(np.amax(coordinates[:, i] + 0.00001))

    coordinates = coordinates *  unit.nanometers

    positions = coordinates

    topology.setUnitCellDimensions(max_coord * unit.nanometers)

    return positions, topology


def _update_ligand_pdb_files(Protein):
    """Updates the ligands pdb files

    Parameters
    ---------------
    Protein

    Returns
    ----------
    Protein
    """

    protein_ligand = mdtraj.load(Protein.pdb_file)

    for ligand in Protein.get_ligand_list():

        try:
            ligand_atoms = protein_ligand.top.select(f'resSeq {ligand.resnum}')
        except Exception:
            ligand_atoms = protein_ligand.top.select(f'resname {ligand.resname}')

        ligand_trj = protein_ligand.atom_slice(ligand_atoms)

        if not ligand.pdb_file:
            
            ligand.pdb_file = f'{ligand.resname}_lgand.pdb'

        ligand_trj.save(ligand.pdb_file, force_overwrite=True)

    return Protein


def parametrize_and_protonate(Protein,
    protein_forcefields,
    water_forcefields,
    ligand_forcefield,
    water_model,
    ph=7.0,
    padding=1.0*unit.nanometers,
    neutralize=False):
    """Parametrize and protonate a Protein and it's ligands with openmm/openmmforcefields
    
    Returns
    ----------
    Protein
    """

    only_protein_pdb = f'{Protein.protein_id}_only_protein.pdb'

    _selections.select_protein_and_ions(Protein.pdb_file, output_file=only_protein_pdb)

    Protein = _update_ligand_pdb_files(Protein)

    ligand_mols = []

    # convert to sdf and protonate at the given ph
    for ligand in Protein.get_ligand_list():
        sdf_file = ligand.pdb_file
        sdf_file = sdf_file.rsplit('.', 1)[0]
        sdf_file += '.sdf'

        _obabel.convert_and_protonate_pdb_to_sdf(ligand.pdb_file,
            sdf_file=sdf_file,
            ph=ph,
            ligand_resname=ligand.resname)

        ligand_mols.append(Molecule.from_file(sdf_file))

    if 'gaff' in ligand_forcefield:
        ligand_ff_generator = GAFFTemplateGenerator(molecules=ligand_mols, forcefield=ligand_forcefield)
    else:
        ligand_ff_generator = SMIRNOFFTemplateGenerator(molecules=ligand_mols, forcefield=ligand_forcefield)

    system_kwargs = {'constraints': None, 'rigidWater': False, 'nonbondedMethod': PME}

    forcefield = ForceField(*protein_forcefields, *water_forcefields)

    forcefield.registerTemplateGenerator(ligand_ff_generator.generator)

    #TODO
    #standard substitution, with attention to CYM

    # Make the protein ligand .top and .pdb file
    protein_ligand_pdb = PDBFile(Protein.pdb_file)
    protein_ligand_modeller = Modeller(protein_ligand_pdb.topology,
        protein_ligand_pdb.positions)

    # Add hydrogens
    protein_ligand_modeller.addHydrogens(forcefield=forcefield, pH=ph, variants=None)

    #TODO custom residue names for Zinc complexing residues
    
    for mol in ligand_mols:
        protein_ligand_modeller.add(mol.to_topology().to_openmm(), mol.conformers[0])

    # Traslate in order to have the minimum of
    # X Y Z to 0 0 0
    protein_ligand_modeller.positions, protein_ligand_modeller.topology =  _scale_positions_and_update_box_vectors(
            protein_ligand_modeller.positions, protein_ligand_modeller.topology)

    protein_ligand_modeller.addSolvent(forcefield,
        model=water_model,
        padding=padding,
        neutralize=neutralize)

    # Write the new protonated PDB
    with open(Protein.pdb_file, 'w') as f:
        PDBFile.writeFile(protein_ligand_modeller.topology, protein_ligand_modeller.positions, file=f)

    # Make the protein ligand .top file
    protein_ligand_system = forcefield.createSystem(
        protein_ligand_modeller.topology, **system_kwargs)

    pmd_structure = parmed.openmm.load_topology(protein_ligand_modeller.topology,
    system=protein_ligand_system, xyz=protein_ligand_modeller.positions)

    Protein.top_file = 'protein_ligand.top'
    pmd_structure.save(Protein.top_file, overwrite=True)

    del protein_ligand_modeller
    del protein_ligand_system

    # Now make .top files for
    # every ligand in vacuum and in water
    # and make a water box on the fly
    water_box = testsystems.WaterBox(box_edge=3.0*unit.nanometers, model=water_model,
                                 constrained=False, nonbondedMethod=PME, ionic_strength=0*unit.molar)

    with open('only_water.pdb', 'w') as f:
        PDBFile.writeFile(water_box.topology, water_box.positions, file=f)

    pmd_structure = parmed.openmm.load_topology(water_box.topology,
    system=water_box.system, xyz=water_box.positions)

    only_water_top = f'only_water_{water_model}.top'
    pmd_structure.save(only_water_top, overwrite=True)

    for lig, mol in zip(Protein.get_ligand_list(), ligand_mols):
        mol_positions, mol_topology =  _scale_positions_and_update_box_vectors(mol.conformers[0],
            mol.to_topology().to_openmm())

        mol_topology.setUnitCellDimensions([5.0, 5.0, 5.0] * unit.nanometers)

        ligand_system = forcefield.createSystem(mol_topology, **system_kwargs)

        pmd_structure = parmed.openmm.load_topology(mol_topology,
        system=ligand_system, xyz=mol_positions)

        lig.top_file = lig.pdb_file.rsplit('.', 1)[0] + '.top'
        pmd_structure.save(lig.top_file, overwrite=True)

        with open(lig.pdb_file, 'w') as f:
            PDBFile.writeFile(mol_topology, mol_positions, file=f)

        # Make the top for water ligand
        # Water first ligand later

        ligand_modeller = Modeller(water_box.topology, water_box.positions)
        ligand_modeller.add(mol_topology, mol_positions)

        ligand_system = forcefield.createSystem(
            ligand_modeller.topology, **system_kwargs)

        pmd_structure = parmed.openmm.load_topology(ligand_modeller.topology,
            system=ligand_system, xyz=ligand_modeller.positions)

        pmd_structure.save('water_' + lig.top_file, overwrite=True)

    #TODO test
    # TODO decide wht to return, a named tuple maybe?
