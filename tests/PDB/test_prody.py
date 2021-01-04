######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################


import unittest
import unittest.mock as mock

from HPC_Drug.PDB import prody


class test_parse_pdb(unittest.TestCase):

    def test_with_correct_input(self):

        with mock.patch("HPC_Drug.PDB.prody.prody.parsePDB", return_value = 'test', autospec=True) as mocked_function:

            output = prody.parse_pdb(file_name = 'file.pdb')

            mocked_function.assert_called_once_with('file.pdb')

            self.assertEqual(output, "test")

class test_write_pdb(unittest.TestCase):

    def test_with_correct_input(self):

        with mock.patch("HPC_Drug.PDB.prody.prody.writePDB", autospec=True) as mocked_function:
                
            prody.write_pdb(structure = "dummy", file_name = 'file.pdb')

            mocked_function.assert_called_once_with('file.pdb', "dummy")


class test_select(unittest.TestCase):

    def test_with_correct_input(self):

        mocked_structure = mock.Mock()

        mocked_structure.select.return_value = "test"

        output = prody.select(structure = mocked_structure, string = "string")

        mocked_structure.select.assert_called_once_with("string")

        self.assertEqual(output, "test")


class test_ProdySelect(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):

        cls.structure = "test"
        
        cls.ProdySelect = prody.ProdySelect(cls.structure)


    def test_only_protein_with_correct_input(self):

        with mock.patch("HPC_Drug.PDB.prody.select", return_value = "test", autospec=True) as mocked_function:

            output = self.ProdySelect.only_protein()

            mocked_function.assert_called_once_with(structure = self.structure, string = "protein")

            self.assertEqual(output, "test")

    def test_protein_and_ions_with_correct_input(self):

        with mock.patch("HPC_Drug.PDB.prody.select", return_value = "test", autospec=True) as mocked_function:

            output = self.ProdySelect.protein_and_ions()

            mocked_function.assert_called_once_with(structure = self.structure, string = "protein or ion")

            self.assertEqual(output, "test")

    def test_resname_with_correct_input(self):

        with mock.patch("HPC_Drug.PDB.prody.select", return_value = "test", autospec=True) as mocked_function:

            output = self.ProdySelect.resname(resname = "resname")

            mocked_function.assert_called_once_with(structure = self.structure, string = "resname RESNAME")

            self.assertEqual(output, "test")

    def test_resnum_with_correct_input(self):

        with mock.patch("HPC_Drug.PDB.prody.select", return_value = "test", autospec=True) as mocked_function:

            output = self.ProdySelect.resnum(resnum = 1)

            mocked_function.assert_called_once_with(structure = self.structure, string = "resnum 1")

            self.assertEqual(output, "test")