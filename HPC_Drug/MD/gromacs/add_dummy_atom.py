######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################
"""Helper function to add a heavy dummy atom
"""

import FSDAMGromacs.gro_files as _gro_files
import FSDAMGromacs.itp_files as _itp_files
import FSDAMGromacs.top_files as _top_files
import PythonPDBStructures.geometry as _geometry
import PythonPDBStructures.pdb.prody_utils as _prody_utils
import PythonPDBStructures.pdb.biopython_utils as _biopython_utils

import HPC_Drug.files_IO.write_on_files as _write_file
import HPC_Drug.auxiliary_functions.path as _path

def add_heavy_dummy_atom(Protein):
    """Adds a heavy dummy atom to a protein gro & top file

    The atom is added on the center of mass of the protein

    Parameters
    -----------
    Protein : HPC_Drug.structures.protein.Protein
    
    Returns
    ------------
    Protein : HPC_Drug.structures.protein.Protein
        the updated protein
    """

    #TODO
    #find protein center of mass
    structure = _prody_utils.parse_pdb(Protein.pdb_file)

    selecter = _prody_utils.ProdySelect(structure)

    structure = selecter.only_protein()

    _prody_utils.write_pdb(
        structure=structure,
        file_name=f'{Protein.protein_id}_only_protein.pdb')

    structure = _biopython_utils.parse_pdb(
        protein_id=Protein.protein_id,
        file_name=f'{Protein.protein_id}_only_protein.pdb'
    )

    COM = _geometry.get_center_of_mass(
        structure=structure,
        geometric=False)

    #Angstrom to nm
    COM = COM * 0.1

    #add the atom to the gro file
    _gro_files.add_atom_to_gro_file(
        input_gro_file=Protein.gro_file,
        output_gro_file=Protein.gro_file,
        coordinates=COM,
        velocities=None,
        atom_name='DU',
        atom_residue_name='DUM'
    )


    #add dummy atom to topology
    atom_types, itp_file = _itp_files.create_dummy_atom_itp(
        mass=1.E20,
        charge=0,
        sigma=0,
        epsilon=0,
        name='DUM',
        charge_group=1,
        atom_number=1
    )

    _write_file.write_file(
        lines=atom_types,
        file_name = "DUM_atomtypes.itp"
    )

    _write_file.write_file(
        lines=itp_file,
        file_name = "DUM.itp"
    )

    atom_types_path = _path.absolute_filepath("DUM_atomtypes.itp")
    itp_file_path = _path.absolute_filepath("DUM.itp")

    _top_files.add_include_after_FF(
        include_line=atom_types_path,
        input_top_file=Protein.top_file,
        output_top_file=Protein.top_file
    )

    _top_files.add_include_after_atomtypes(
        include_line=itp_file_path,
        input_top_file=Protein.top_file,
        output_top_file=Protein.top_file
    )

    _top_files.add_molecules(
        name="DUM",
        number=1,
        input_top_file=Protein.top_file,
        output_top_file=Protein.top_file
    )

    return Protein