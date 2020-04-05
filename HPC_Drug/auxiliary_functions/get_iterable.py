######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains the function to get an iterable
"""

import collections


def get_iterable(x):
    """Returns an iterable, even if given a single value
    if x is a string returns (string,) even though a string is an iterable"""
    
    if type(x) == str:

        return (x,)

    elif isinstance(x, collections.Iterable):

        return x

    else:

        return (x,)