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

from HPC_Drug.PDB import mmcif2pdb

class test_mmcif2pdb(unittest.TestCase):

    def test_with_TypeError(self):

        with self.assertRaises(TypeError):

            Protein = mock.Mock()
            Protein.file_type = "WRONG_TYPE"

            mmcif2pdb.mmcif2pdb(Protein = Protein)

    def test_with_pdb(self):

        Protein = mock.Mock()
        Protein.file_type = "pdb"
        Protein.pdb_file = "test.pdb"

        Protein = mmcif2pdb.mmcif2pdb(Protein = Protein)

        Protein.update_structure.assert_not_called()

        self.assertEqual(Protein.pdb_file, "test.pdb")
        self.assertEqual(Protein.file_type, "pdb")

    def test_with_cif(self):

        with mock.patch("HPC_Drug.PDB.biopython.write_mmcif") as mocked_mmcfi2pdb:

            Protein = mock.Mock()
            Protein.file_type = "cif"
            Protein.pdb_file = "test.cif"

            Protein = mmcif2pdb.mmcif2pdb(Protein = Protein)

            Protein.update_structure.assert_called_once()

            mocked_mmcfi2pdb.assert_called_once()

            self.assertEqual(Protein.pdb_file, "test.pdb")
            self.assertEqual(Protein.file_type, "pdb")