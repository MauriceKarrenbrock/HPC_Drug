######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains the generic functions needed to write or append lines to a non formatted text file (.txt .pdb .tex ...)
"""

def write_file(lines, file_name = "file.txt"):
    """This function writes a new file (or overwrides an existing one)
    
    lines :: anything that supports iteration like lists, tuples etc or a single string

    file_name :: string, default "new_file.txt"
    
    doesn't return anything
    
    The string must be correctly formatted with all the needed newlines when given as input"""

    if type(lines) == str:
        lines = [lines]

    with open(file_name, 'w') as f:

        for line in lines:
            
            f.write(line)

def append_file(lines, file_name = "file.txt"):
    """This function appendsa lines to an existing file
    
    lines :: anything that supports iteration like lists, tuples etc or a single string

    file_name :: string, default "new_file.txt"
    
    doesn't return anything
    
    The string must be correctly formatted with all the needed newlines when given as input"""

    if type(lines) == str:
        lines = [lines]

    with open(file_name, 'a') as f:

        for line in lines:
            
            f.write(line)