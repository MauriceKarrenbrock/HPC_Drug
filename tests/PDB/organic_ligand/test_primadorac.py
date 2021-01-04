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

from HPC_Drug.PDB.organic_ligand import primadorac

class test_primadorac(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        Protein = mock.Mock()

        lig_1 = mock.Mock()
        #lig_2 = mock.Mock()

        ligands = [lig_1]

        for i in range(len(ligands)):

            ligands[i].resname = "LIG"
            ligands[i].pdb_file = "pdb_file"
            ligands[i].itp_file = "itp_file"
            ligands[i].tpg_file = "tpg_file"
            ligands[i].prm_file = "prm_file"


        Protein.get_ligand_list.return_value = ligands

        with mock.patch("HPC_Drug.auxiliary_functions.path.absolute_programpath", return_value = "primadorac_path"):

            cls.primadorac = primadorac.Primadorac(Protein = Protein, primadorac_path = "primadorac_path", ph = 7.0)



    def test__init__(self):

        with mock.patch("HPC_Drug.auxiliary_functions.path.absolute_programpath", return_value = "test_path") as mocked_abspath:

            output = primadorac.Primadorac(Protein = "test_protein", primadorac_path = "test_path", ph = 0)

            mocked_abspath.assert_called_once_with(program = "test_path")

            self.assertEqual(output.Protein, "test_protein")
            self.assertEqual(output.primadorac_path, "test_path")
            self.assertEqual(output.ph, 0)


    def test__run(self):

        with mock.patch("HPC_Drug.auxiliary_functions.run.subprocess_run") as mocked_run:

            self.primadorac._run(string = "string")

            mocked_run.assert_called_once_with(commands = "string",
                                            shell = False,
                                            error_string = "primadorac failure")


    def test__rename_itp_file_exists(self):

        with mock.patch("HPC_Drug.PDB.organic_ligand.primadorac.os.path.exists", return_value = True) as mocked_exists:

            output = self.primadorac._rename_itp(file_to_search = "test_file", ligand_resname = "LIG")

            mocked_exists.assert_called_once()

            self.assertEqual(output, "test_file")

    
    def test__rename_itp_file_itp(self):
        
        with mock.patch("HPC_Drug.PDB.organic_ligand.primadorac.os.path.exists", side_effect = [False,True]) as mocked_exists:

            with mock.patch("HPC_Drug.PDB.organic_ligand.primadorac.os.rename", return_value = "test.itp") as mocked_rename:

                output = self.primadorac._rename_itp(file_to_search = "test_file", ligand_resname = "LIG")

                mocked_exists.assert_called()

                mocked_rename.assert_called_once_with("file.itp", "LIG.itp")

                self.assertEqual(output, "LIG.itp")


    def test__rename_itp_file_itp_None(self):
        
        with mock.patch("HPC_Drug.PDB.organic_ligand.primadorac.os.path.exists", side_effect = [False,False,True]) as mocked_exists:

            with mock.patch("HPC_Drug.PDB.organic_ligand.primadorac.os.rename", return_value = "test.itp") as mocked_rename:

                output = self.primadorac._rename_itp(file_to_search = "test_file", ligand_resname = "LIG")

                mocked_exists.assert_called()

                mocked_rename.assert_called_once_with('file.itp none', 'LIG.itp')

                self.assertEqual(output, "LIG.itp")


    def test__rename_itp_RuntimeError(self):

        with mock.patch("HPC_Drug.PDB.organic_ligand.primadorac.os.path.exists", side_effect = [False,False,False]) as mocked_exists:

            with self.assertRaises(RuntimeError):

                self.primadorac._rename_itp(file_to_search = "test_file", ligand_resname = "LIG")
                
                mocked_exists.assert_called()


    def test__edit_itp(self):

        lines = [
            "1\n",
            "2\n",
            "3\n",
            "4\n",
            "5",
            "6",
            "7",
            "8",
            "9",
            "\n\n\naaa LIG aaaa\n\n\n",
            "\n\nbbsbsb bsbbs name-p dvsd\n"
        ]

        expected_output_lines = [
            "aaa TEST_RESNAME aaaa\n",
            "bbsbsb bsbbs TEST_RESNAME dvsd\n"
        ]

        with mock.patch("HPC_Drug.files_IO.read_file.read_file", return_value = lines) as mocked_read:

            with mock.patch("HPC_Drug.files_IO.write_on_files.write_file") as mocked_write:

                output = self.primadorac._edit_itp(ligand_resname = "TEST_RESNAME", itp_file = "test_file")

                mocked_read.assert_called_once_with(file_name = "test_file")

                mocked_write.assert_called_once_with(lines = expected_output_lines, file_name = "test_file")

                self.assertEqual(output, "test_file")

    
    def test_execute(self):

        def dummy(path):
            return path

        with mock.patch("HPC_Drug.auxiliary_functions.path.absolute_programpath", return_value = "bash") as mocked_program_path:

            with mock.patch("HPC_Drug.auxiliary_functions.path.absolute_filepath", side_effect = dummy) as mocked_file_path:

                with mock.patch.object(primadorac.Primadorac, "_run") as mocked_run:

                    with mock.patch.object(primadorac.Primadorac, "_edit_itp") as mocked_edit:

                        with mock.patch.object(primadorac.Primadorac, "_rename_itp") as mocked_rename:

                            output = self.primadorac.execute()

                            mocked_program_path.assert_called_once_with(program = "bash")

                            mocked_file_path.assert_called()

                            mocked_run.assert_called_with(string = [
                                "bash",
                                self.primadorac.primadorac_path,
                                '-gp',
                                "pdb_file"])

                            mocked_edit.assert_called()

                            mocked_rename.assert_called_with(
                                file_to_search = 'pdb_file' + '-p' + '.itp',
                                ligand_resname = self.primadorac.Protein.get_ligand_list()[-1].resname
                            )


            


