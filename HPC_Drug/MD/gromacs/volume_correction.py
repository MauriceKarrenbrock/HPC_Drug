# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""function to calculate the volume correction from *_pullx.xvg files
"""

from pathlib import Path
import numpy as np

from FSDAMGromacs.parse import GromacsParsePullDistances

from PythonFSDAM.free_energy_calculations import volume_correction as _volume_correction

def volume_correction(directory=None, temperature= 298.15):
    """Calculates the volume correction from the various BATTERY* dirs

    Will look for ny *_pullx.xvg file in any BATTARY*/scaled0/ directory
    and use all the distances taken from there to calculate the correction
    in Kcal/mol
    
    Parameters
    -------------
    directory : Path or str, default current directory
        the directory in which the BATTERY* dirs are
    temperature : float, default=298.15
        Kelvin

    Returns
    ---------
    float
        free energy correction in Kcal/mol
        shall be added to the dissociation free energy
        or subtracted from the binding free energy

    Raises
    -----------
    RuntimeError
        if no BATTERY*/scaled0/*_pullx.xvg is found in the bound dir
    """

    distances = np.array([])

    if directory is None:
        directory = '.'
    directory = Path(directory)

    # Should be only one file per BATTERY but I want to make
    # this function as generic as possible
    files = list(directory.glob('BATTERY*/scaled0/*_pullx.xvg'))

    if not files:
        raise RuntimeError("No BATTERY*/scaled0/*_pullx.xvg found in the bound dir")

    for ff in files:
        distances_dict = GromacsParsePullDistances.parse(str(ff.resolve()))

        distances = np.concatenate((distances, distances_dict[1]))

        # Free memory (it is a big matrix)
        distances_dict = None

    return _volume_correction(distance_values=distances, temperature=temperature)