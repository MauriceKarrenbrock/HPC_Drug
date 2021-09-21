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

import warnings

import numpy as np

from simtk import unit
from simtk.openmm.app import PDBFile, ForceField, PME, Modeller
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator, SMIRNOFFTemplateGenerator
from openmmtools import testsystems

import parmed
import mdtraj

import HPC_Drug.PDB.complex_selections as _selections
import HPC_Drugs.MD.openbabel_utils as _obabel
from HPC_Drug.MD import residue_renaming

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
    neutralize=False,
    residue_renaming_type='standard'):
    """Parametrize and protonate a Protein and it's ligands with openmm/openmmforcefields

    Parameters
    --------------
    Protein : HPC_Drug.structures.protein.Protein
    protein_forcefields : list(str)
        the strings must be the forcefields names that simtk.openmm.Forcefield can find
        therefore something in openmm, openmmforcefields or in the working directory
    water_forcefields : list(str)
        the strings must be the forcefields names that simtk.openmm.Forcefield can find
        therefore something in openmm, openmmforcefields or in the working directory
    ligand_forcefield : str
        one of the FF supported by openmmforcefields.generators, at the moment GAFFx.x and the
        openff initiative ones (es SMIRNOFF)
    water_model : str
        a water model like tip3p, it must be supported by openmm.Modeller (usuallu only the
        classic 3 point ones)
    ph : float, default=7.0
        at which ph to protonate the ligand
    padding : simtik.unit.Quantity or float, default=1.0*unit.nanometers
        how much water padding to do, if a float is given it will be considered in nm
    neutralize : bool, default=False
    residue_renaming_type str or bool, default=standard
        the kind of residue renaming for the given FF
        if False no residue will be renamed
    
    Returns
    ----------
    Protein, only_water_pdb, only_water_top
        each ligand of the Protein has an updated top_file and solvated_top_file
    """

    if not unit.is_quantity(padding):
        padding = padding * unit.nanometers

    # I first do standard because the protonation part
    # does not recognize custom stuff
    if residue_renaming_type:
        if 'amber' in (''.join(protein_forcefields)).lower():
            Protein = residue_renaming.ResidueRenamer(
                        Protein=Protein,
                        forcefield='amber',
                        substitution="standard",
                        ph=ph).execute()
        else:
            warnings.warn('only amber FF residue renamings are supported at the moment, therefore '
            'carefully check the output!')

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

    # Make the protein ligand .top and .pdb file
    protein_ligand_pdb = PDBFile(only_protein_pdb)
    protein_ligand_modeller = Modeller(protein_ligand_pdb.topology,
        protein_ligand_pdb.positions)

    for mol in ligand_mols:
        protein_ligand_modeller.add(mol.to_topology().to_openmm(), mol.conformers[0])


    # Create variants
    variants = []
    tmp_top = mdtraj.load(Protein.pdb_file).topology
    for residue in tmp_top.residues:
        if residue.name.upper() == 'CYM':
            variants.append('CYX')
        else:
            variants.append(None)


    # Add hydrogens
    protein_ligand_modeller.addHydrogens(forcefield=forcefield, pH=ph, variants=variants)

    # Custom residue names for example for Zinc complexing residues
    # but to do it I will have to redo a lot of stuff
    # TODO refactor it
    if residue_renaming_type != 'standard' and residue_renaming_type:

        # Temporarely write the pdb file to be modified
        with open(Protein.pdb_file, 'w') as f:
            PDBFile.writeFile(protein_ligand_modeller.topology, protein_ligand_modeller.positions, file=f)

        if 'amber' in (''.join(protein_forcefields)).lower():
            Protein = residue_renaming.ResidueRenamer(
                    Protein=Protein,
                    forcefield='amber',
                    substitution=residue_renaming_type,
                    ph=ph).execute()

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

            # Make the protein ligand .top and .pdb file
            protein_ligand_pdb = PDBFile(only_protein_pdb)
            protein_ligand_modeller = Modeller(protein_ligand_pdb.topology,
                protein_ligand_pdb.positions)

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

    only_water_pdb = f'only_water_{water_model}.pdb'
    with open(only_water_pdb, 'w') as f:
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

        lig.solvated_top_file = 'water_' + lig.top_file
        pmd_structure.save(lig.solvated_top_file, overwrite=True)

    return Protein, only_water_pdb, only_water_top
