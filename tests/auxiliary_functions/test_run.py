######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

import unittest
import unittest.mock as mock

from HPC_Drug.auxiliary_functions import run

class test_subrocess_run(unittest.TestCase):

    def test_with_returncode_0(self):

        mocked_returncode = mock.Mock()
        mocked_returncode.returncode = 0
        mocked_returncode.stdout = ""
        mocked_returncode.stderr = ""

        with mock.patch("HPC_Drug.auxiliary_functions.run.subprocess.run", return_value = mocked_returncode) as mocked_subprocess_run:

            with mock.patch("HPC_Drug.auxiliary_functions.run.subprocess.PIPE") as mocked_pipe:

                with mock.patch("HPC_Drug.auxiliary_functions.run.print") as mocked_print:

                    run.subprocess_run(commands = "test")

                    mocked_subprocess_run.assert_called_once_with("test",
                                                                shell = False,
                                                                stdout=mocked_pipe,
                                                                stderr=mocked_pipe,
                                                                universal_newlines = False)

                    mocked_print.assert_called()

        
    def test_with_returncode_not_0(self):

        mocked_returncode = mock.Mock()
        mocked_returncode.returncode = -1
        mocked_returncode.stdout = ""
        mocked_returncode.stderr = ""

        with mock.patch("HPC_Drug.auxiliary_functions.run.subprocess.run", return_value = mocked_returncode) as mocked_subprocess_run:

            with mock.patch("HPC_Drug.auxiliary_functions.run.subprocess.PIPE") as mocked_pipe:

                with mock.patch("HPC_Drug.auxiliary_functions.run.print") as mocked_print:

                    with self.assertRaises(RuntimeError):
                    
                        run.subprocess_run(commands = "test")

                        mocked_subprocess_run.assert_called_once_with("test",
                                                                    shell = False,
                                                                    stdout=mocked_pipe,
                                                                    stderr=mocked_pipe,
                                                                    universal_newlines = False)

                        mocked_print.assert_called()

