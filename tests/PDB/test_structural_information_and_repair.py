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

import HPC_Drug.PDB.structural_information_and_repair as struct_rep

class test_InfoRepair(unittest.TestCase):

    def test__init__(self):

        output = struct_rep.InfoRepair(Protein = "protein", repairing_method = "pdbfixer")

        self.assertEqual(output.Protein, "protein")
        self.assertEqual(output.repairing_method, "pdbfixer")

    def test__parse_header(self):

        with mock.patch("HPC_Drug.PDB.structural_information.mmcif_header.get_metalbinding_disulf_ligands", return_value = ("protein", "list")):

            output = struct_rep.InfoRepair(Protein = "test", repairing_method = "pdbfixer")

            output._parse_header()

            self.assertEqual(output.Protein, "protein")
            self.assertEqual(output.organic_ligand_list, "list")

    def test__parse_structure(self):

        with mock.patch("HPC_Drug.PDB.structural_information.scan_structure.get_metalbinding_disulf_ligands", return_value = ("protein", "list")):

            output = struct_rep.InfoRepair(Protein = "test", repairing_method = "pdbfixer")

            output._parse_structure()

            self.assertEqual(output.Protein, "protein")
            self.assertEqual(output.organic_ligand_list, "list")

    def test__repair(self):

        with mock.patch("HPC_Drug.PDB.repair_pdb.repair.repair", return_value = "repaired", autospec = True):

            output = struct_rep.InfoRepair(Protein = "test", repairing_method = "pdbfixer")

            output._repair()

            self.assertEqual(output.Protein, "repaired")

    def test__pdb(self):

        with mock.patch.object(struct_rep.InfoRepair, "_repair", autospec = True) as mocked_repair:

            with mock.patch.object(struct_rep.InfoRepair, "_parse_structure", autospec = True) as mocked_parse:

                output = struct_rep.InfoRepair(Protein = "test", repairing_method = "pdbfixer")

                output._pdb()

                mocked_repair.assert_called_once()
                mocked_parse.assert_called_once()

    def test_cif_header(self):

        with mock.patch.object(struct_rep.InfoRepair, "_repair", autospec = True) as mocked_repair:

            with mock.patch.object(struct_rep.InfoRepair, "_parse_header", autospec = True) as mocked_parse:

                output = struct_rep.InfoRepair(Protein = "test", repairing_method = "pdbfixer")

                output._cif()

                mocked_repair.assert_called_once()
                mocked_parse.assert_called_once()

    def test_cif_structure(self):

        with mock.patch.object(struct_rep.InfoRepair, "_repair", autospec = True) as mocked_repair:

            with mock.patch.object(struct_rep.InfoRepair, "_parse_header", side_effect = Exception, autospec = True) as mocked_parse:

                with mock.patch.object(struct_rep.InfoRepair, "_pdb", autospec = True) as mocked_pdb:

                    output = struct_rep.InfoRepair(Protein = "test", repairing_method = "pdbfixer")

                    output._cif()

                    mocked_repair.assert_not_called()
                    mocked_parse.assert_called_once()
                    mocked_pdb.assert_called_once()

    def test_get_info_and_repair_typeerror(self):

        Protein = mock.Mock()
        Protein.file_type = "WRONG_TYPE"
        
        output = struct_rep.InfoRepair(Protein = Protein, repairing_method = "pdbfixer")

        with self.assertRaises(TypeError):

            output.get_info_and_repair()

    def test_get_info_and_repair_cif(self):

        with mock.patch.object(struct_rep.InfoRepair, "_cif", return_value = ("output", "list"), autospec = True) as mocked_cif:

            Protein = mock.Mock()
            Protein.file_type = "cif"
            
            output = struct_rep.InfoRepair(Protein = Protein, repairing_method = "pdbfixer")

            out_out = output.get_info_and_repair()

            mocked_cif.assert_called_once()

            self.assertEqual(out_out, (output.Protein, output.organic_ligand_list))

    def test_get_info_and_repair_pdb(self):

        with mock.patch.object(struct_rep.InfoRepair, "_pdb", return_value = ("output", "list"), autospec = True) as mocked_pdb:

            Protein = mock.Mock()
            Protein.file_type = "pdb"
            
            output = struct_rep.InfoRepair(Protein = Protein, repairing_method = "pdbfixer")

            out_out = output.get_info_and_repair()

            mocked_pdb.assert_called_once()

            self.assertEqual(out_out, (output.Protein, output.organic_ligand_list))