"""
Copyright (c) 2020-2020 Maurice Karrenbrock

This software is open-source and is distributed under the
GNU Affero General Public License v3 (agpl v3) license

A copy of the license must be included with any copy of the program or part of it
"""

import unittest
import unittest.mock as mock

from HPC_Drug.PDB.repair_pdb import repair

class test_repair(unittest.TestCase):

    def test_with_non_existing_repairing_method(self):

        with self.assertRaises(NotImplementedError):

            repair.repair(Protein = "test", repairing_method = "WRONG_METHOD")

    def test_with_correct_input(self):

        with mock.patch("HPC_Drug.PDB.repair_pdb.pdbfixer.repair", return_value = "repaired") as mocked_repair:

            output = repair.repair(Protein = "test", repairing_method = "pdbfixer")

            mocked_repair.assert_called_once_with(Protein = "test")

            self.assertEqual(output, "repaired")
