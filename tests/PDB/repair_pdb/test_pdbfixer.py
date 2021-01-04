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

from HPC_Drug.PDB.repair_pdb import pdbfixer

class test__repair(unittest.TestCase):
    """Tests _repair the private function"""

    def test_with_wrong_file_type(self):

        with self.assertRaises(TypeError):

            pdbfixer._repair(input_file_name = "test_input",
                            file_type = "WRONG",
                            output_file_name = "test_output")

    def test_with_add_H_False_pdb(self):

        with mock.patch("HPC_Drug.PDB.repair_pdb.pdbfixer.open") as mocked_context_manager:

            with mock.patch("HPC_Drug.PDB.repair_pdb.pdbfixer.pdbfixer.PDBFixer", autospec = True) as mocked_PDBFixer:

                with mock.patch("HPC_Drug.PDB.repair_pdb.pdbfixer.simtk.openmm.app.PDBFile.writeFile") as mocked_pdbwriter:
                
                    mocked_context_manager.return_value.__enter__.return_value = "StringIO"

                    mocked_fixer = mock.Mock()

                    mocked_PDBFixer.return_value = mocked_fixer


                    #The function to test
                    output = pdbfixer._repair(input_file_name = "test_input",
                                    file_type = "pdb",
                                    output_file_name = "test_output",
                                    add_H = False)

                    
                    mocked_context_manager.assert_called()

                    mocked_fixer.findMissingResidues.assert_called_once()

                    mocked_fixer.findNonstandardResidues.assert_called_once()

                    mocked_fixer.replaceNonstandardResidues.assert_called_once()

                    mocked_fixer.findMissingAtoms.assert_called_once()

                    mocked_fixer.addMissingAtoms.assert_called_once()

                    self.assertEqual(output, "test_output")


    def test_with_add_H_True_pdb(self):

        with mock.patch("HPC_Drug.PDB.repair_pdb.pdbfixer.open") as mocked_context_manager:

            with mock.patch("HPC_Drug.PDB.repair_pdb.pdbfixer.pdbfixer.PDBFixer", autospec = True) as mocked_PDBFixer:

                with mock.patch("HPC_Drug.PDB.repair_pdb.pdbfixer.simtk.openmm.app.PDBFile.writeFile") as mocked_pdbwriter:
                
                    mocked_context_manager.return_value.__enter__.return_value = "StringIO"

                    mocked_fixer = mock.Mock()

                    mocked_PDBFixer.return_value = mocked_fixer


                    #The function to test
                    output = pdbfixer._repair(input_file_name = "test_input",
                                    file_type = "pdb",
                                    output_file_name = "test_output",
                                    add_H = True,
                                    ph = 7.0)

                    
                    mocked_context_manager.assert_called()

                    mocked_fixer.findMissingResidues.assert_called_once()

                    mocked_fixer.findNonstandardResidues.assert_called_once()

                    mocked_fixer.replaceNonstandardResidues.assert_called_once()

                    mocked_fixer.findMissingAtoms.assert_called_once()

                    mocked_fixer.addMissingAtoms.assert_called_once()

                    mocked_fixer.addMissingHydrogens.assert_called_once_with(7.0)

                    mocked_pdbwriter.assert_called()

                    self.assertEqual(output, "test_output")

    #WON'T BE TESTED AS LONG AS THERE IS THE patch
    # def test_with_add_H_True_cif(self):

    #     with mock.patch("HPC_Drug.PDB.repair_pdb.pdbfixer.open") as mocked_context_manager:

    #         with mock.patch("HPC_Drug.PDB.repair_pdb.pdbfixer.pdbfixer.PDBFixer", autospec = True) as mocked_PDBFixer:

    #             with mock.patch("HPC_Drug.PDB.repair_pdb.pdbfixer.simtk.openmm.app.pdbxfile.PDBxFile.writeFile", autospec = True) as mocked_pdbwriter:

    #                 mocked_context_manager.return_value.__enter__.return_value = "StringIO"

    #                 mocked_fixer = mock.Mock()

    #                 mocked_PDBFixer.return_value = mocked_fixer
            
    #                 #The function to test
    #                 output = pdbfixer._repair(input_file_name = "test_input",
    #                                 file_type = "cif",
    #                                 output_file_name = "test_output",
    #                                 add_H = True,
    #                                 ph = 7.0)

                    
    #                 mocked_context_manager.assert_called()

    #                 mocked_fixer.findMissingResidues.assert_called_once()

    #                 mocked_fixer.findNonstandardResidues.assert_called_once()

    #                 mocked_fixer.replaceNonstandardResidues.assert_called_once()

    #                 mocked_fixer.findMissingAtoms.assert_called_once()

    #                 mocked_fixer.addMissingAtoms.assert_called_once()

    #                 mocked_fixer.addMissingHydrogens.assert_called_once_with(7.0)

    #                 mocked_pdbwriter.assert_called()

    #                 self.assertEqual(output, "test_output")


class test_repair(unittest.TestCase):

    def test_with_correct_input(self):

        def dummy(path):
            return path

        with mock.patch("HPC_Drug.auxiliary_functions.path.absolute_filepath", side_effect = dummy):

            with mock.patch("HPC_Drug.PDB.repair_pdb.pdbfixer._repair", return_value = "file") as mocked_repair:

                Protein = mock.Mock()
                Protein.pdb_file = None

                Protein = pdbfixer.repair(Protein = Protein)

                mocked_repair.assert_called_once()

                self.assertEqual(Protein.pdb_file, "file")