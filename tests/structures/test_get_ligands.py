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

from HPC_Drug.structures import get_ligands

class test_get_ligands(unittest.TestCase):

    def test_reise_typeError(self):

        Protein = mock.Mock()

        Protein.file_type = "WRONG_TYPE"

        with self.assertRaises(TypeError):

            get_ligands.get_ligands(Protein = Protein, ligand_resnames_resnums = [["LIG", 1]])

    def test_empty_liglist(self):

        Protein = mock.Mock()

        Protein.file_type = "pdb"

        output = get_ligands.get_ligands(Protein = Protein, ligand_resnames_resnums = [])

        output.add_ligand.assert_not_called()

        self.assertEqual(Protein, output)

    def test_None_liglist(self):

        Protein = mock.Mock()

        Protein.file_type = "pdb"

        output = get_ligands.get_ligands(Protein = Protein, ligand_resnames_resnums = None)

        output.add_ligand.assert_not_called()

        self.assertEqual(Protein, output)

    def test_with_valid_liglist(self):
        #should be done better

        with mock.patch("HPC_Drug.structures.ligand.Ligand", autospec = True) as mocked_lig:

            with mock.patch("HPC_Drug.PDB.prody.ProdySelect", autospec = True) as mocked_select:

                mocked_lig.return_value = mocked_lig

                mocked_lig.resnum = 0
                mocked_lig.resname = "LIG"

                Protein = mock.Mock()

                Protein.file_type = "pdb"

                get_ligands.get_ligands(Protein = Protein, ligand_resnames_resnums = [["LIG", 1], ["LIG", 2]])

                mocked_lig.assert_called()

                mocked_select.assert_called()

                Protein.add_ligand.assert_called()