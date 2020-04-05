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

from HPC_Drug.PDB import biopython


class test_parse_pdb(unittest.TestCase):

    def test_with_correct_input(self):
        with mock.patch('HPC_Drug.PDB.biopython.Bio.PDB.PDBParser', autospec=True) as mocked_class:

            mocked_class.return_value.get_structure.return_value = 'test_structure'

            output = biopython.parse_pdb(protein_id = '2f3z', file_name = 'tests/files4tests/2f3z.pdb')

            mocked_class.assert_called_once()

            mocked_class.return_value.get_structure.assert_called_once_with('2f3z', 'tests/files4tests/2f3z.pdb')

            self.assertEqual('test_structure', output)

    def test_with_wrong_input_type(self):

        wrong_inputs = ((1, 'tests/files4tests/2f3z.pdb'), ('2f3z', 7), (1.1, 'tests/files4tests/2f3z.pdb'), ('2f3z', 7.7))

        with self.assertRaises(TypeError):

            for i in wrong_inputs:

                biopython.parse_pdb(protein_id = i[0], file_name = i[1])



class test_parse_mmcif(unittest.TestCase):

    def test_with_correct_input(self):
        with mock.patch('HPC_Drug.PDB.biopython.Bio.PDB.MMCIFParser', autospec=True) as mocked_class:

            mocked_class.return_value.get_structure.return_value = 'test_structure'

            output = biopython.parse_mmcif(protein_id = '2f3z', file_name = 'tests/files4tests/2f3z.cif')

            mocked_class.assert_called_once()

            mocked_class.return_value.get_structure.assert_called_once_with('2f3z', 'tests/files4tests/2f3z.cif')

            self.assertEqual('test_structure', output)

    def test_with_wrong_input_type(self):

        wrong_inputs = ((1, 'tests/files4tests/2f3z.cif'), ('2f3z', 7), (1.1, 'tests/files4tests/2f3z.cif'), ('2f3z', 7.7))

        with self.assertRaises(TypeError):

            for i in wrong_inputs:

                biopython.parse_mmcif(protein_id = i[0], file_name = i[1])


class test_parse_structure_factory(unittest.TestCase):

    def test_with_wrong_input(self):

        with mock.patch('HPC_Drug.structures.protein.Protein', file_type = "WRONG", autospec=True) as mocked_protein:

            with self.assertRaises(TypeError):

                biopython.structure_factory(Protein = mocked_protein)

    def test_with_pdb(self):

        with mock.patch('HPC_Drug.structures.protein.Protein', file_type = "pdb", pdb_file = "file.pdb", protein_id = "idid", autospec=True) as mocked_protein:

            with mock.patch('HPC_Drug.PDB.biopython.parse_pdb', return_value = "test") as mocked_parser:

                output = biopython.structure_factory(mocked_protein)

                mocked_parser.assert_called_once_with(protein_id = 'idid', file_name = 'file.pdb')

                self.assertEqual(output, "test")

    def test_with_mmcif(self):

        with mock.patch('HPC_Drug.structures.protein.Protein', file_type = "cif", pdb_file = "file.cif", protein_id = "idid", autospec=True) as mocked_protein:

            with mock.patch('HPC_Drug.PDB.biopython.parse_mmcif', return_value = "test") as mocked_parser:

                output = biopython.structure_factory(mocked_protein)

                mocked_parser.assert_called_once_with(protein_id = 'idid', file_name = 'file.cif')

                self.assertEqual(output, "test")

        
class test_parse_mmcif2dict(unittest.TestCase):

    def test_with_wrong_input(self):

        wrong_inputs = (1, 1.1, list())
        with self.assertRaises(TypeError):

            for i in wrong_inputs:
                biopython.mmcif2dict(i)

    def test_with_correct_input(self):

        with mock.patch("HPC_Drug.PDB.biopython.Bio.PDB.MMCIF2Dict.MMCIF2Dict", return_value = dict(), autospec=True) as mocked_function:

            output = biopython.mmcif2dict("file.cif")

            mocked_function.assert_called_once_with("file.cif")

            self.assertEqual(output, dict())


class test_write_pdb(unittest.TestCase):

    def test_with_correct_input(self):

        with mock.patch('HPC_Drug.PDB.biopython.isinstance', return_value = True):

            with mock.patch('HPC_Drug.PDB.biopython.Bio.PDB.PDBIO', autospec=True) as mocked_class:
                
                MockedStructure = mock.Mock()

                biopython.write_pdb(structure = MockedStructure, file_name = 'file.pdb')

                mocked_class.assert_called_once()

                mocked_class.return_value.set_structure.assert_called_once_with(MockedStructure)

                mocked_class.return_value.save.assert_called_once_with('file.pdb')


    def test_with_wrong_structure_type(self):

        with self.assertRaises(TypeError):

            wrong_inputs = (1, 1.1, "aaa", list(), tuple(), dict())

            for i in wrong_inputs:
            
                biopython.write_pdb(structure = i, file_name = "file.pdb")

    def test_with_wrong_file_name_type(self):

        with self.assertRaises(TypeError):

            with mock.patch('HPC_Drug.PDB.biopython.isinstance', return_value = True):

                wrong_inputs = (1, 1.1, list(), tuple(), dict())

                for i in wrong_inputs:
                
                    biopython.write_pdb(structure = "dummy", file_name = i)



class test_write_mmcif(unittest.TestCase):

    def test_with_correct_input(self):

        with mock.patch('HPC_Drug.PDB.biopython.isinstance', return_value = True):

            with mock.patch('HPC_Drug.PDB.biopython.Bio.PDB.MMCIFIO', autospec=True) as mocked_class:
                
                MockedStructure = mock.Mock()

                biopython.write_mmcif(structure = MockedStructure, file_name = 'file.cif')

                mocked_class.assert_called_once()

                mocked_class.return_value.set_structure.assert_called_once_with(MockedStructure)

                mocked_class.return_value.save.assert_called_once_with('file.cif')


    def test_with_wrong_structure_type(self):

        with self.assertRaises(TypeError):

            wrong_inputs = (1, 1.1, "aaa", list(), tuple(), dict())

            for i in wrong_inputs:
            
                biopython.write_mmcif(structure = i, file_name = "file.cif")

    def test_with_wrong_file_name_type(self):

        with self.assertRaises(TypeError):

            with mock.patch('HPC_Drug.PDB.biopython.isinstance', return_value = True):

                wrong_inputs = (1, 1.1, list(), tuple(), dict())

                for i in wrong_inputs:
                
                    biopython.write_mmcif(structure = "dummy", file_name = i)


class test_write_dict2mmcif(unittest.TestCase):

    def test_with_correct_input(self):

        with mock.patch('HPC_Drug.PDB.biopython.Bio.PDB.MMCIFIO', autospec=True) as mocked_class:

            MockedDictionary = mock.Mock()

            biopython.write_dict2mmcif(dictionary = MockedDictionary, file_name = 'file.cif')

            mocked_class.assert_called_once()

            mocked_class.return_value.set_dict.assert_called_once_with(MockedDictionary)

            mocked_class.return_value.save.assert_called_once_with('file.cif')

    def test_with_wrong_file_name_type(self):

        with self.assertRaises(TypeError):

            wrong_inputs = (1, 1.1, list(), tuple(), dict())

            for i in wrong_inputs:
            
                biopython.write_dict2mmcif(dictionary = "dummy", file_name = i)

    # def test_with_wrong_dictionary_type(self):

    #     with self.assertRaises(TypeError):

    #         wrong_inputs = (1, 1.1)

    #         for i in wrong_inputs:
            
    #             biopython.write_dict2mmcif(dictionary = i, file_name = "tests/files4tests/file.cif")

    #     with self.assertRaises(UnboundLocalError):
            
    #         wrong_inputs = (list(), tuple())

    #         for i in wrong_inputs:

    #             biopython.write_dict2mmcif(dictionary = i, file_name = "tests/files4tests/file.cif")
        
        
    #     with self.assertRaises(ValueError):

    #         biopython.write_dict2mmcif(dictionary = "aaa", file_name = "tests/files4tests/file.cif")


class test_write(unittest.TestCase):

    def test_with_wrong_input(self):

        with self.assertRaises(TypeError):

            biopython.write(structure = "dummy", file_type = "WRONG", file_name = None)


    def test_with_dict_and_pdb(self):

        with self.assertRaises(RuntimeError):

            biopython.write(structure = dict(), file_type = "pdb", file_name = None)

    def test_with_dict_and_mmcif(self):

        with mock.patch('HPC_Drug.PDB.biopython.write_dict2mmcif') as mocked_function:

            biopython.write(structure = dict(), file_type = "cif", file_name = "file.cif")

            mocked_function.assert_called_once_with(dictionary = dict(), file_name = "file.cif")


    def test_with_structure_and_mmcif(self):

        with mock.patch('HPC_Drug.PDB.biopython.write_mmcif') as mocked_function:

            biopython.write(structure = "dummy", file_type = "cif", file_name = "file.cif")

            mocked_function.assert_called_once_with(structure = "dummy", file_name = "file.cif")

    def test_with_structure_and_pdb(self):

        with mock.patch('HPC_Drug.PDB.biopython.write_pdb') as mocked_function:

            biopython.write(structure = "dummy", file_type = "pdb", file_name = "file.pdb")

            mocked_function.assert_called_once_with(structure = "dummy", file_name = "file.pdb")