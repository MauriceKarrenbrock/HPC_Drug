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

from HPC_Drug.PDB import merge_pdb

class test_merge_pdb(unittest.TestCase):

    def test_with_correct_input(self):

        input_list = [
            ["ATOM   1286  CA  LYS B   1      20.249  -7.881   8.457  1.00 65.74           C\n",
            "ATOM   1286  CA  LYS B   2      20.249  -7.881   8.457  1.00 65.74           C\n",
            "CONNETC1286  CA  LYS B  -1      20.249  -7.881   8.457  1.00 65.74           C\n"],

            ["ATOM   1286  CA  LYS B   1      20.249  -7.881   8.457  1.00 65.74           C\n",
            "HETATM 1286  CA  LYS B   2      20.249  -7.881   8.457  1.00 65.74           C\n",
            "CONNETC1286  CA  LYS B  -1      20.249  -7.881   8.457  1.00 65.74           C\n"],

            ["ATOM   1286  CA  LYS B   1      20.249  -7.881   8.457  1.00 65.74           C\n",
            "TER    1286  CA  LYS B   2      20.249  -7.881   8.457  1.00 65.74           C\n",
            "CONNETC1286  CA  LYS B  -1      20.249  -7.881   8.457  1.00 65.74           C\n"]
        ]

        expected_output = [
            ["ATOM   1286  CA  LYS B   1      20.249  -7.881   8.457  1.00 65.74           C\n",
            "ATOM   1286  CA  LYS B   2      20.249  -7.881   8.457  1.00 65.74           C\n",
            "HETATM   86  CA  LIG B   3      20.000  -7.881   8.457  1.00 65.74           C\n",
            "HETATM   86  CA  LIG B   4      20.000  -7.881   8.457  1.00 65.74           C\n",
            "CONNETC1286  CA  LYS B  -1      20.249  -7.881   8.457  1.00 65.74           C\n"],

            ["ATOM   1286  CA  LYS B   1      20.249  -7.881   8.457  1.00 65.74           C\n",
            "HETATM 1286  CA  LYS B   2      20.249  -7.881   8.457  1.00 65.74           C\n",
            "HETATM   86  CA  LIG B   3      20.000  -7.881   8.457  1.00 65.74           C\n",
            "HETATM   86  CA  LIG B   4      20.000  -7.881   8.457  1.00 65.74           C\n",
            "CONNETC1286  CA  LYS B  -1      20.249  -7.881   8.457  1.00 65.74           C\n"],

            ["ATOM   1286  CA  LYS B   1      20.249  -7.881   8.457  1.00 65.74           C\n",
            "TER    1286  CA  LYS B   2      20.249  -7.881   8.457  1.00 65.74           C\n",
            "HETATM   86  CA  LIG B   3      20.000  -7.881   8.457  1.00 65.74           C\n",
            "HETATM   86  CA  LIG B   4      20.000  -7.881   8.457  1.00 65.74           C\n",
            "CONNETC1286  CA  LYS B  -1      20.249  -7.881   8.457  1.00 65.74           C\n"]
        ]

        for i in range(len(input_list)):

            with mock.patch("HPC_Drug.files_IO.read_file.read_file", side_effect = [input_list[i],\
                ["HETATM   86  CA  LIG B0000      20.000  -7.881   8.457  1.00 65.74           C\n"],\
                    ["HETATM   86  CA  LIG B0000      20.000  -7.881   8.457  1.00 65.74           C\n"]]) as mocked_read:

                with mock.patch("HPC_Drug.files_IO.write_on_files.write_file") as mocked_write:

                    with mock.patch("HPC_Drug.auxiliary_functions.path.absolute_filepath", return_value = "abs_path") as mocked_path:

                        Protein = mock.Mock()
                        ligand_1 = mock.Mock()
                        ligand_2 = mock.Mock()

                        ligand_1.pdb_file = "ligand_file"
                        ligand_2.pdb_file = "ligand_file"

                        ligand_1.resnum = 0
                        ligand_2.resnum = 0

                        Protein.protein_id = "test_id"
                        Protein.pdb_file = "protein_file"
                        Protein.get_ligand_list.return_value = [ligand_1, ligand_2]


                        output = merge_pdb.merge_pdb(Protein = Protein)

                        mocked_read.assert_called()

                        mocked_write.assert_called_once_with(lines = expected_output[i], file_name = "test_id_joined.pdb")

                        mocked_path.assert_called_once_with(path = "test_id_joined.pdb")

                        self.assertEqual(output.get_ligand_list()[0].resnum, 3)
                        self.assertEqual(output.get_ligand_list()[1].resnum, 4)
                        self.assertEqual(output.pdb_file, "abs_path")
