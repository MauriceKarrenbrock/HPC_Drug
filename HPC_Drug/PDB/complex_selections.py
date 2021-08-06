######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################
"""Functions to do selection operations that can't be done with a mdtraj selection string"""

import mdtraj
from mdtraj.core.residue_names import _SOLVENT_TYPES, _WATER_RESIDUES

import HPC_Drug.important_lists as _important_lists

def select_protein_and_ions(input_file, output_file, **kwargs):
    """Select protein and all ions from a given file

    Parameters
    -----------
    input_file : str or path
        any file type supported by mdtraj
    output_file : str or path
        any file type supported by mdtraj
    kwargs
        any extra keyword argument that might be needed
        for the mdtraj.load function
    """

    traj = mdtraj.load(str(input_file), **kwargs)

    protein_atoms = list(traj.topology.select('protein'))

    ions_residues = _SOLVENT_TYPES - _WATER_RESIDUES

    ions_atoms = [atom.index for atom in traj.topology.atoms if
                atom.residue.name in ions_residues]

    traj.atom_slice(protein_atoms + ions_atoms, inplace=True)

    traj.save(str(output_file), force_overwrite=True)


def remove_trash(input_file, output_file, trash=None, **kwargs):
    """Removes unwanted residues from a structure given the unwanted residue names

    Parameters
    -----------
    input_file : str or path
        any file type supported by mdtraj
    output_file : str or path
        any file type supported by mdtraj
    trash : iterable(str), optional
        the residue names to remove
        as default it is list(HPC_Drug.important_lists.trash) + list(HPC_Drug.important_lists.trash_ions)
    kwargs
        any extra keyword argument that might be needed
        for the mdtraj.load function
    """

    if trash is None:
        trash = list(_important_lists.trash) + list(_important_lists.trash_ions)

    traj = mdtraj.load(str(input_file), **kwargs)

    good_atoms = [atom.index for atom in traj.topology.atoms if
                atom.residue.name not in trash if
                atom.residue.name.upper() not in trash]

    traj.atom_slice(good_atoms, inplace=True)

    traj.save(str(output_file), force_overwrite=True)