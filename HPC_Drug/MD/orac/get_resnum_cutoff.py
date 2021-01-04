######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
A function used by other functions
"""

def get_first_resnum(structure):
    """
    returns the resnumber of the first residue

    structure :: biopython structure
    """

    r = structure.get_residues()

    for i in r:
        resnum = i.id[1]
        break
    
    return resnum

def get_resnum_cutoff(structure):
    """
    returns 1 - resnum of the first residue

    useful with Orac that starts counting residues by one

    structure :: biopython structure
    """

    cutoff = 1 - get_first_resnum(structure = structure)

    return cutoff