######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################


import unittest
import unittest.mock

import HPC_Drug.pipelines

class test_Pipeline(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        # make a dummy __init__() function
        with unittest.mock.patch.object(HPC_Drug.pipelines.Pipeline(), "__init__", lambda x: None):

            cls.pipeline = HPC_Drug.pipelines.Pipeline()


    def test_get_protein_file_exception_wrong_local_input(self):

        with self.assertRaises(ValueError):

            self.pipeline.local = "WRONG"

            self.pipeline.get_protein_file()

    def _dummy_path(self, path):
        return path

    
    def test_get_protein_file_with_local_yes_existing_file(self):

        with unittest.mock.patch("HPC_Drug.pipelines.os.path.exists", return_value = True):

            with unittest.mock.patch("HPC_Drug.auxiliary_functions.path.absolute_filepath", side_effect = self._dummy_path):
            
                self.pipeline.local = 'yes'

                self.pipeline.protein_filename = 'test_name'

                self.pipeline.get_protein_file()

                self.assertEqual('test_name', self.pipeline.protein_filename)

    def test_get_protein_file_with_local_yes_not_existing_file(self):
        
            with unittest.mock.patch("HPC_Drug.pipelines.os.path.exists", return_value = False):

                with self.assertRaises(FileNotFoundError):

                    self.pipeline.local = 'yes'

                    self.pipeline.get_protein_file()

    def test_get_protein_file_with_local_no(self):

        with unittest.mock.patch("HPC_Drug.PDB.download_pdb.download", return_value = "test", autospec=True) as mocked_function:

            with unittest.mock.patch("HPC_Drug.auxiliary_functions.path.absolute_filepath", side_effect = self._dummy_path):
            
                self.pipeline.local = 'no'

                self.pipeline.protein_filetype = "cif"

                self.pipeline.protein_id = "idid"

                self.pipeline.get_protein_file()

                mocked_function.assert_called_once_with(protein_id = self.pipeline.protein_id, file_type = self.pipeline.protein_filetype, pdir = None)

                self.assertEqual(self.pipeline.protein_filename, 'test')