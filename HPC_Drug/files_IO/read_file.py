######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains the function to read a standard .txt file
"""

def read_file(file_name):
    """Reads the given text file and returns a list
    containing the lines
    
    file :: str or path
    
    Can have problems with very big files"""

    with open(file_name, 'r') as f:

        lines = f.readlines()

    return lines

