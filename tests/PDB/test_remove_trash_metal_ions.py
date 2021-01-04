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

from HPC_Drug.PDB import remove_trash_metal_ions


class test_remove_trash_metal_ions(unittest.TestCase):

    #define a dummy trash list for testing purposes
    trash = ["TRS"]

    def test_raise_TypeError(self):

        Protein = mock.Mock()

        Protein.file_type = "WRONG_TYPE"

        Protein.pdb_file = "test"

        with mock.patch("HPC_Drug.files_IO.read_file.read_file") as mocked_read:

            with self.assertRaises(TypeError):

                remove_trash_metal_ions.remove_trash_metal_ions(Protein = Protein, trash = self.trash)

                mocked_read.assert_called_once_with(file_name = "test")

    
    def test_pdb(self):

        Protein = mock.Mock()

        Protein.file_type = "pdb"

        Protein.pdb_file = "test"

        lines = [
            "HETATM 1461  O   TRS A 331      -3.404   5.794  21.117  1.00 35.49           O  ",
            "ATOM   1462  O   TRS A 332      -3.404   5.794  21.117  1.00 35.49           O  ",
            "TER    1463  O   TRS A 333      -3.404   5.794  21.117  1.00 35.49           O  ",
            "HETATM 1464  O   HOH A 334      -3.404   5.794  21.117  1.00 35.49           O  ",
            "ATOM   1465  O   CYS A 335      -3.404   5.794  21.117  1.00 35.49           O  ",
            "TER    1466  O   HIS A 336      -3.404   5.794  21.117  1.00 35.49           O  "
        ]

        expected_output_list = [
            "HETATM 1464  O   HOH A 334      -3.404   5.794  21.117  1.00 35.49           O\n",
            "ATOM   1465  O   CYS A 335      -3.404   5.794  21.117  1.00 35.49           O\n",
            "TER    1466  O   HIS A 336      -3.404   5.794  21.117  1.00 35.49           O\n",
            "end\n"
        ]

        with mock.patch("HPC_Drug.files_IO.read_file.read_file", return_value = lines) as mocked_read:

            with mock.patch("HPC_Drug.files_IO.write_on_files.write_file") as mocked_write:

                output = remove_trash_metal_ions.remove_trash_metal_ions(Protein = Protein, trash = self.trash)

                mocked_read.assert_called_once_with(file_name = "test")

                mocked_write.assert_called_once_with(lines = expected_output_list, file_name = "test")

                self.assertEqual(output.pdb_file, "test")


    def test_cif(self):

        Protein = mock.Mock()

        Protein.file_type = "cif"

        Protein.pdb_file = "test"

        lines = [
            "# ",
            "loop_",
            "_atom_site.group_PDB ",
            "_atom_site.id ",
            "_atom_site.type_symbol ",
            "_atom_site.label_atom_id ",
            "_atom_site.label_alt_id ",
            "ATOM   1    N  N   . TRS A 1 3   ? 12.469  7.601   -5.979 1.00 51.86 ? 3    TRS A N   1",
            "HETATM 2    N  N   . TRS A 1 3   ? 12.469  7.601   -5.979 1.00 51.86 ? 3    TRS A N   1",
            "ATOM   3    N  N   . GLN A 1 3   ? 12.469  7.601   -5.979 1.00 51.86 ? 3    GLN A N   1",
            "HETATM 4    N  N   . GLN A 1 3   ? 12.469  7.601   -5.979 1.00 51.86 ? 3    GLN A N   1"
        ]

        expected_output_list = [
            "#\n",
            "loop_\n",
            "_atom_site.group_PDB\n",
            "_atom_site.id\n",
            "_atom_site.type_symbol\n",
            "_atom_site.label_atom_id\n",
            "_atom_site.label_alt_id\n",
            "ATOM   3    N  N   . GLN A 1 3   ? 12.469  7.601   -5.979 1.00 51.86 ? 3    GLN A N   1\n",
            "HETATM 4    N  N   . GLN A 1 3   ? 12.469  7.601   -5.979 1.00 51.86 ? 3    GLN A N   1\n"
        ]

        with mock.patch("HPC_Drug.files_IO.read_file.read_file", return_value = lines) as mocked_read:

            with mock.patch("HPC_Drug.files_IO.write_on_files.write_file") as mocked_write:

                output = remove_trash_metal_ions.remove_trash_metal_ions(Protein = Protein, trash = self.trash)

                mocked_read.assert_called_once_with(file_name = "test")

                mocked_write.assert_called_once_with(lines = expected_output_list, file_name = "test")

                self.assertEqual(output.pdb_file, "test")
