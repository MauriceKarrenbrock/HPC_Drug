######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################


import unittest

from HPC_Drug.files_IO import read_file

class test_write_file(unittest.TestCase):

    def test_wrong_input(self):

        wrong_inputs = (1, 1.5)

        with self.assertRaises(TypeError):

            for i in wrong_inputs:
                
                read_file.read_file(file_name = i)

    def test_with_file(self):

        file_name = "tests/files4tests/file_read_file.txt"

        lines = read_file.read_file(file_name = file_name)

        self.assertEqual(lines, ["aaa\n", "BBB\n", "111\n", "2.2.2\n"])
