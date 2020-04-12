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

from HPC_Drug.MD.orac import orac_input


class test_OracInput(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        Protein = mock.Mock()
        ligand_1 = mock.Mock()
        ligand_2 = mock.Mock()

        ligand_1.resname = "LIG"
        ligand_2.resname = "LIG"

        ligand_1.resnum = 1
        ligand_2.resnum = 2

        ligand_1.pdb_file = "ligand_pdb"
        ligand_2.pdb_file = "ligand_pdb"

        ligand_1.structure = "lig_struct"
        ligand_2.structure = "lig_struct"

        ligand_1.prm_file = "ligand_prm"
        ligand_2.prm_file = "ligand_prm"

        ligand_1.tpg_file = "ligand_tpg"
        ligand_2.tpg_file = "ligand_tpg"

        Protein.protein_id = "id"
        Protein.pdb_file = "protein_file"
        Protein.prm_file = "prm_file"
        Protein.tpg_file = "tpg_file"
        Protein.structure = "protein_structure"
        Protein.sulf_bonds = [(1,2),(3,4)]

        Protein.get_ligand_list.return_value = [ligand_1, ligand_2]

        with mock.patch.object(orac_input.OracInput, "__init__", lambda x, Protein: None):

            cls.orac_object = orac_input.OracInput( Protein = "dummy")

            cls.orac_object.Protein = Protein

            cls.orac_object.MD_program_path = "orac"

            cls.orac_object.solvent_pdb = "water.pdb"

            cls.orac_object.orac_in_file = "id_orac.in"

            cls.orac_object.output_pdb_file = "id_orac.pdb"

            cls.orac_object.template = []

            orient = mock.Mock()

            cls.orac_object.orient = orient


    def test__init__(self):

        with mock.patch("HPC_Drug.orient.Orient", return_value = "test_orient", autospec = True) as mocked_orient:

            with mock.patch("HPC_Drug.auxiliary_functions.path.absolute_programpath", return_value = "abs_orac", autospec = True) as mocked_programpath:

                with mock.patch("HPC_Drug.MD.orac.orac_input.os.getcwd", return_value = "work_dir", autospec = True) as mocked_workdir:

                    with mock.patch("HPC_Drug.MD.orac.orac_input.importlib_resources.path", autospec = True) as mocked_importlib:

                        with mock.patch("HPC_Drug.MD.orac.orac_input.str", return_value = "water.pdb") as mocked_str_cast:


                            Protein = mock.Mock()

                            Protein.protein_id = "id"

                            test_class = orac_input.OracInput(Protein = Protein, solvent_pdb = None, MD_program_path = "orac")

                            mocked_orient.assert_called_once()

                            mocked_programpath.assert_called_once_with(program = "orac")

                            mocked_workdir.assert_called()

                            mocked_importlib.assert_any_call("HPC_Drug.lib", "water.pdb")

                            mocked_str_cast.assert_called_once()

                            self.assertEqual(test_class.Protein, Protein)
                            self.assertEqual(test_class.orac_in_file, "work_dir" + "/id" + "_orac.in")
                            self.assertEqual(test_class.solvent_pdb, "water.pdb")
                            self.assertEqual(test_class.MD_program_path, "abs_orac")
                            self.assertEqual(test_class.output_pdb_file, "work_dir" + "/id" + "_orac.pdb")
                            self.assertEqual(test_class.template, [])
                            self.assertEqual(test_class.orient, "test_orient")


    def test_write_box(self):

        self.orac_object.orient.create_box.return_value = (1.000000, 2.00000000, 3.00000000)

        output = self.orac_object._write_box()

        self.assertEqual(output, "   CRYSTAL    1.00    2.00    3.00  !! simulation BOX")

    def test__write_sulf_bond_string(self):

        with mock.patch.object(orac_input.OracInput, "_get_protein_resnumber_cutoff", return_value = 0, autospec = True) as mocked_cutoff:

            output = self.orac_object._write_sulf_bond_string()

            mocked_cutoff.assert_called()

            self.assertEqual(output, "   bond 1sg 2sg residue   1  2\n   bond 1sg 2sg residue   2  1\n   bond 1sg 2sg residue   3  4\n   bond 1sg 2sg residue   4  3")


    def test__get_protein_resnumber_cutoff(self):

        residue = mock.Mock()
        residue.id = ("", 5, "")

        structure = mock.Mock()
        structure.get_residues.return_value = [residue]

        self.orac_object.Protein.structure = structure

        output = self.orac_object._get_protein_resnumber_cutoff()

        self.assertEqual(output, -4)

        #reset to standard value
        self.orac_object.Protein.structure = "protein_structure"


    def test__write_ligand_tpg_path(self):
        
        output = self.orac_object._write_ligand_tpg_path()

        self.assertEqual(output, "   READ_TPG_ASCII ligand_tpg !! ligand\n   READ_TPG_ASCII ligand_tpg !! ligand")

    
    def test__write_ligand_prm_path(self):
        
        output = self.orac_object._write_ligand_prm_path()

        self.assertEqual(output, "   READ_PRM_ASCII ligand_prm  !! ligand\n   READ_PRM_ASCII ligand_prm  !! ligand")


    def test__get_ligand_name_from_tpg(self):

        lines = [
            "RANDOM AAAAA BBBBB ZZZZZZ",
            "RESIDUE LIG ZZZZZZZZZZZ",
            "11111111111111111"
        ]

        with mock.patch("HPC_Drug.files_IO.read_file.read_file", return_value = lines) as mocked_read:
        
            output = self.orac_object._get_ligand_name_from_tpg()

            mocked_read.assert_called()

            self.assertEqual(output, "      LIG !! ligand name in tpg file\n      LIG !! ligand name in tpg file")


    def test__get_ligand_name_from_tpg_RuntimeError(self):

        lines = [
            "RANDOM AAAAA BBBBB ZZZZZZ",
            "222222222 ZZZZZZZZZZZ",
            "11111111111111111"
        ]

        with mock.patch("HPC_Drug.files_IO.read_file.read_file", return_value = lines):

            with self.assertRaises(RuntimeError):
        
                self.orac_object._get_ligand_name_from_tpg()


    def test__write_solvent_grid(self):

        self.orac_object.orient.reset_mock()

        self.orac_object.orient.create_box.return_value = (1, 2, 3)

        self.orac_object.orient.create_solvent_grid.return_value = (1, 2, 3)

        output = self.orac_object._write_solvent_grid()

        self.orac_object.orient.create_box.assert_called_once()

        self.orac_object.orient.create_solvent_grid.assert_called_once()

        self.assertEqual(output, "   GENERATE   1  2  3  !! the grid depends on the BOX")


    def test__write_EWALD_PME(self):

        self.orac_object.orient.reset_mock()

        self.orac_object.orient.create_box.return_value = (1, 2, 3)

        self.orac_object.orient.create_recipr_solvent_grid.return_value = (1, 2, 3)

        output = self.orac_object._write_EWALD_PME()

        self.orac_object.orient.create_box.assert_called_once()

        self.orac_object.orient.create_recipr_solvent_grid.assert_called_once()

        self.assertEqual(output, "   EWALD PME 0.37   1  2  3   4 !! grid on the reciprocal")


    #def test__write_ADD_STR_COM(self):
