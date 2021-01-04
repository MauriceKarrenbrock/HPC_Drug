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
import os

from HPC_Drug.PDB import download_pdb


class test_download(unittest.TestCase):

    def test_with_correct_input_pdb(self):

        with mock.patch("HPC_Drug.PDB.download_pdb.os.path.exists", return_value = True):

            with mock.patch("HPC_Drug.PDB.download_pdb.Bio.PDB.PDBList", autospec=True) as mocked_class:

                mocked_class.return_value.retrieve_pdb_file.return_value = "test"

                output = download_pdb.download(protein_id = "idid", file_type = "pdb", pdir = None)

                mocked_class.assert_called_once_with()

                mocked_class.return_value.retrieve_pdb_file.assert_called_once_with("idid", False, os.getcwd(), file_format = "pdb", overwrite = True)

                self.assertEqual(output, "test")

    def test_with_correct_input_cif(self):

        with mock.patch("HPC_Drug.PDB.download_pdb.os.path.exists", return_value = True):

            with mock.patch("HPC_Drug.PDB.download_pdb.Bio.PDB.PDBList", autospec=True) as mocked_class:

                mocked_class.return_value.retrieve_pdb_file.return_value = "test"

                for i in ("cif", "mmCif"):

                    mocked_class.reset_mock()

                    output = download_pdb.download(protein_id = "idid", file_type = i, pdir = None)

                    mocked_class.assert_called_once_with()

                    mocked_class.return_value.retrieve_pdb_file.assert_called_once_with("idid", False, os.getcwd(), file_format = "mmCif", overwrite = True)

                    self.assertEqual(output, "test")

    def test_download_did_not_happen(self):

        with mock.patch("HPC_Drug.PDB.download_pdb.os.path.exists", return_value = False):

            with mock.patch("HPC_Drug.PDB.download_pdb.Bio.PDB.PDBList", autospec=True):

                with self.assertRaises(FileNotFoundError):
                
                    download_pdb.download(protein_id = "idid", file_type = "pdb", pdir = None)