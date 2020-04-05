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

import HPC_Drug.structures.structure as structure

class test_Structure(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):

        with mock.patch.object(structure.Structure, "__init__", return_value = None):

            cls.structure = structure.Structure()

    def test_init_raise(self):

        with self.assertRaises(NotImplementedError):
            
            structure.Structure()

    def test_write_cif_prody(self):

        self.structure.file_type = "cif"

        self.structure.pdb_file = "test_file"

        self.structure.structure = "test_structure"

        with self.assertRaises(TypeError):

            self.structure.write(struct_type = 'prody')

    def _dummy_path(self, path):
        return path

    def test_write_pdb_prody(self):

        self.structure.file_type = "pdb"

        self.structure.pdb_file = "test_file"

        self.structure.structure = "test_structure"

        with mock.patch("HPC_Drug.PDB.prody.write_pdb", autospec = True) as mocked_function:

            with mock.patch("HPC_Drug.auxiliary_functions.path.absolute_filepath", side_effect = self._dummy_path):

                self.structure.write(struct_type = 'prody')

                mocked_function.assert_called_once_with(structure = self.structure.structure, file_name = self.structure.pdb_file)

    def test_write_biopython(self):

        self.structure.file_type = "pdb"

        self.structure.pdb_file = "test_file"

        self.structure.structure = "test_structure"

        with mock.patch("HPC_Drug.PDB.biopython.write", autospec = True) as mocked_function:

            with mock.patch("HPC_Drug.auxiliary_functions.path.absolute_filepath", side_effect = self._dummy_path):

                self.structure.write(struct_type = 'biopython')

                mocked_function.assert_called_once_with(structure = self.structure.structure, file_type = self.structure.file_type, file_name = self.structure.pdb_file)

    

    def test_update_structure_prody_cif(self):

        self.structure.file_type = "cif"

        with self.assertRaises(TypeError):

            self.structure.update_structure(struct_type = "prody")


    def test_update_structure_prody_pdb(self):

        self.structure.file_type = "pdb"

        self.structure.pdb_file = "test_file"

        self.structure.structure = "test_structure"

        with mock.patch("HPC_Drug.PDB.prody.parse_pdb", return_value = "NEW_test_structure", autospec = True) as mocked_function:

            self.structure.update_structure(struct_type = "prody")

            mocked_function.assert_called_once_with(file_name = self.structure.pdb_file)

            self.assertEqual(self.structure.structure, "NEW_test_structure")