######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################
"""Helper function to add a heavy dummy atom
"""

import mdtraj

import FSDAMGromacs.gro_files as _gro_files
import FSDAMGromacs.itp_files as _itp_files
import FSDAMGromacs.top_files as _top_files


def add_heavy_dummy_atom(gro_file,
    top_file,
    output_gro='with_dummy.gro',
    output_top='with_dummy.top'):
    """Adds a heavy dummy atom to a protein gro & top file

    The atom is added on the center of mass of the protein

    Parameters
    -----------
    gro_file : str
    top_file : str
    output_gro : str, default='with_dummy.gro'
    output_top : str, default='with_dummy.top'
    
    Returns
    ------------
    output_gro, output_top : str, str
    """

    #find protein center of mass in nm
    traj = mdtraj.load(gro_file)

    COM = mdtraj.geometry.distance.compute_center_of_mass(traj, select='protein')[0]

    #add the atom to the gro file
    _gro_files.add_atom_to_gro_file(
        input_gro_file=gro_file,
        output_gro_file=output_gro,
        coordinates=COM,
        velocities=(0., 0., 0.),
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

    _top_files.add_string_after_FF(
        input_string=''.join(atom_types),
        input_top_file=top_file,
        output_top_file=output_top
    )

    _top_files.add_string_after_atomtypes(
        input_string=''.join(itp_file),
        input_top_file=output_top,
        output_top_file=output_top
    )

    _top_files.add_molecules(
        name="DUM",
        number=1,
        input_top_file=output_top,
        output_top_file=output_top
    )

    return output_gro, output_top