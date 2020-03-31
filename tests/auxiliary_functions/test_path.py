"""
Copyright (c) 2020-2020 Maurice Karrenbrock

This software is open-source and is distributed under the
GNU Affero General Public License v3 (agpl v3) license

A copy of the license must be included with any copy of the program or part of it
"""

import unittest
import unittest.mock as mock

from HPC_Drug.auxiliary_functions import path



class test_absolute_filepath(unittest.TestCase):

    def test_file_not_found(self):

        with mock.patch("HPC_Drug.auxiliary_functions.path.os.path.exists", return_value = False, autospec = True):

            with self.assertRaises(FileNotFoundError):

                path.absolute_filepath(path = "NOT_EXIST")

    def test_existing_file(self):

        with mock.patch("HPC_Drug.auxiliary_functions.path.os.path.exists", return_value = True, autospec = True):

            with mock.patch("HPC_Drug.auxiliary_functions.path.os.path.abspath", return_value = "abs_path", autospec = True) as mocked_abspath:

                with mock.patch("HPC_Drug.auxiliary_functions.path.os.path.expanduser", autospec = True) as mocked_expanduser:

                    with mock.patch("HPC_Drug.auxiliary_functions.path.os.path.expandvars", autospec = True) as mocked_expandvars:
            
                        output = path.absolute_filepath(path = "test")

                        mocked_abspath.assert_called_once()

                        mocked_expanduser.assert_called_once()

                        mocked_expandvars.assert_called_once_with("test")

                        self.assertEqual(output, "abs_path")


class test_which(unittest.TestCase):

    def test_with_non_existing_executable(self):

        with mock.patch("HPC_Drug.auxiliary_functions.path.shutil.which", return_value = None, autospec = True):

            with self.assertRaises(OSError):

                path.which(program = "NOT_EXIST")

    def test_with_existing_executable(self):

        with mock.patch("HPC_Drug.auxiliary_functions.path.shutil.which", return_value = "abs_path", autospec = True) as mocked_function:

            output = path.which(program = "test")

            mocked_function.assert_called_once_with("test")

            self.assertEqual(output, "abs_path")


class test_absolute_programpath(unittest.TestCase):

    def test_raise_OSError(self):

        with mock.patch("HPC_Drug.auxiliary_functions.path.which", side_effect = OSError, autospec = True):

            with mock.patch("HPC_Drug.auxiliary_functions.path.absolute_filepath", side_effect = FileNotFoundError, autospec = True):

                with self.assertRaises(OSError):

                    path.absolute_programpath(program = "test")


    def test_use_absolute_filepath(self):

        with mock.patch("HPC_Drug.auxiliary_functions.path.which", side_effect = OSError, autospec = True):

            with mock.patch("HPC_Drug.auxiliary_functions.path.absolute_filepath", return_value = "test_output", autospec = True) as mocked_function:
                
                output = path.absolute_programpath(program = "test")

                mocked_function.assert_called_once_with(path = "test")

                self.assertEqual(output, "test_output")

    def test_use_which(self):

        with mock.patch("HPC_Drug.auxiliary_functions.path.which", return_value = "test_output", autospec = True) as mocked_function:

            output = path.absolute_programpath(program = "test")

            mocked_function.assert_called_once_with(program = "test")

            self.assertEqual(output, "test_output")