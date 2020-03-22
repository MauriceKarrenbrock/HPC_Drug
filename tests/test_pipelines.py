"""
Copyright (c) 2020-2020 Maurice Karrenbrock

This software is open-source and is distributed under the
GNU Affero General Public License v3 (agpl v3) license

A copy of the license must be included with any copy of the program or part of it
"""

import unittest
import unittest.mock

import HPC_Drug.pipelines

class test_Pipeline(unittest.TestCase):

    def test_download_exception_wrong_local_input(self):

        wrong_names = ('ggg', 'K', 'sissi')

        # make a dummy __init__() function
        with unittest.mock.patch.object(HPC_Drug.pipelines.Pipeline(), "__init__", lambda x: None):
            p = HPC_Drug.pipelines.Pipeline()
            
            for name in wrong_names:
                with self.assertRaises(ValueError):

                    p.local = name

                    name = p.download()

    
    def test_download_with_local_yes_existing_file(self):

        # make a dummy __init__() function
        with unittest.mock.patch.object(HPC_Drug.pipelines.Pipeline(), "__init__", lambda x: None):
            with unittest.mock.patch("os.path.exists", lambda x: True):

                p = HPC_Drug.pipelines.Pipeline()
                
                p.local = 'yes'

                p.protein_filename = 'dummy_name'

                self.assertEqual('dummy_name', p.download())

    def test_download_with_local_yes_not_existing_file(self):

        # make a dummy __init__() function
        with unittest.mock.patch.object(HPC_Drug.pipelines.Pipeline(), "__init__", lambda x: None):
            with unittest.mock.patch("os.path.exists", lambda x: False):
                with unittest.mock.patch("HPC_Drug.file_manipulation.download_protein_structure", lambda x, y, z: 'hello'):

                    p = HPC_Drug.pipelines.Pipeline()
                    
                    p.local = 'yes'

                    self.assertEqual('hello', p.download())

                    



    def test_download_with_local_no(self):

        # make a dummy __init__() function
        with unittest.mock.patch.object(HPC_Drug.pipelines.Pipeline(), "__init__", lambda x: None):
            with unittest.mock.patch("HPC_Drug.file_manipulation.download_protein_structure", lambda x, y, z: 'hello'):
                
                p = HPC_Drug.pipelines.Pipeline()

                p.local = 'no'

                self.assertEqual('hello', p.download())