######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains the functions needed to get absolute paths for files and executables
"""

import os.path
import shutil

def absolute_filepath(path):
    """Takes a string and returns the absolute path of the file
    If the file does not exist raises a FileNotFoundError

    path :: string

    return absolute_path
    """

    absolute_path = os.path.abspath(os.path.expanduser(os.path.expandvars(path)))

    if not os.path.exists(absolute_path):
        raise FileNotFoundError(f"{absolute_path}")

    return absolute_path

def which(program):
    """Uses shutil.which to get the absolute path of an executable that is in path
    example: path = which("python"), path = /usr/bin/python
    If the executable doesn't exist raises a OSError


    program :: string

    If you don't know if the program is in $PATH or not use absolute_programpath(program)
    """

    absolute_path = shutil.which(program)

    if absolute_path == None:
        raise OSError(f"{program} is not in $PATH")

    else:
        return absolute_path


def absolute_programpath(program):
    """It returns the absolute path to a program both if it is in $PATH or not

    if the executable doesn't exist raises an OSError

    program :: string
    """

    try:

        try:

            absolute_path = which(program = program)

        except OSError:

            absolute_path = absolute_filepath(path = program)

    except FileNotFoundError:

        raise OSError(f"{program} doesn't exist both in PATH and in the normal directories")

    return absolute_path