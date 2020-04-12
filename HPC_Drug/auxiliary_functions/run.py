######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains the functions needed to run external programs
"""

import subprocess
import os

def subprocess_run(commands,
    shell = False,
    universal_newlines = False,
    error_string = "error during the call of an external program",
    cwd = os.getcwd()):

    """
    runs an external program using subprocess.run

    if it fails will print the standard output, standard error
    and raise RuntimeError

    commands :: list , it is the list of strings containing the command that
    subprocess.run will run

    shell :: bool, dafault False, if == True the commands will be executed in the shell

    universa_newlines :: bool, dafault False

    error_string :: the sting to give to the RuntimeError as argument

    cwd :: string , default the current working directory, it is the working directory for the child process
    """

    r = subprocess.run(commands,
                    shell = shell,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines = universal_newlines,
                    cwd = cwd)

    print(r.stdout)
    print(r.stderr)

    if r.returncode != 0:
        raise RuntimeError(error_string)