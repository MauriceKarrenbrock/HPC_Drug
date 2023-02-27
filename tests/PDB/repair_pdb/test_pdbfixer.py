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

class test_repair(unittest.TestCase):
    """Tests repair the private function"""

    def test_with_wrong_file_type(self):

        with self.assertRaises(TypeError):

            pdbfixer.repair(input_file_name = "test_input.WRONG",
                            output_file_name = "test_output")

    def test_with_add_H_False_pdb(self):

        with mock.patch("HPC_Drug.PDB.repair_pdb.pdbfixer.open") as mocked_context_manager:

            with mock.patch("HPC_Drug.PDB.repair_pdb.pdbfixer.pdbfixer.PDBFixer", autospec = True) as mocked_PDBFixer:

                with mock.patch("HPC_Drug.PDB.repair_pdb.pdbfixer.simtk.openmm.app.PDBFile.writeFile") as mocked_pdbwriter:
                
                    mocked_context_manager.return_value.__enter__.return_value = "StringIO"

                    mocked_fixer = mock.Mock()

                    mocked_PDBFixer.return_value = mocked_fixer


                    #The function to test
                    output = pdbfixer.repair(input_file_name = "test_input.pdb",
                                    output_file_name = "test_output.pdb",
                                    add_H = False)

                    
                    mocked_context_manager.assert_called()

                    mocked_fixer.findMissingResidues.assert_called_once()

                    mocked_fixer.findNonstandardResidues.assert_called_once()

                    mocked_fixer.replaceNonstandardResidues.assert_called_once()

                    mocked_fixer.findMissingAtoms.assert_called_once()

                    mocked_fixer.addMissingAtoms.assert_called_once()

                    mocked_pdbwriter.assert_called()

                    self.assertEqual(output, "test_output.pdb")


    def test_with_add_H_True_pdb(self):

        with mock.patch("HPC_Drug.PDB.repair_pdb.pdbfixer.open") as mocked_context_manager:

            with mock.patch("HPC_Drug.PDB.repair_pdb.pdbfixer.pdbfixer.PDBFixer", autospec = True) as mocked_PDBFixer:

                with mock.patch("HPC_Drug.PDB.repair_pdb.pdbfixer.simtk.openmm.app.PDBFile.writeFile") as mocked_pdbwriter:
                
                    mocked_context_manager.return_value.__enter__.return_value = "StringIO"

                    mocked_fixer = mock.Mock()

                    mocked_PDBFixer.return_value = mocked_fixer


                    #The function to test
                    output = pdbfixer.repair(input_file_name = "test_input.pdb",
                                    output_file_name = "test_output.pdb",
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

                    self.assertEqual(output, "test_output.pdb")

    #WON'T BE TESTED AS LONG AS THERE IS THE patch
    # def test_with_add_H_True_cif(self):

    #     with mock.patch("HPC_Drug.PDB.repair_pdb.pdbfixer.open") as mocked_context_manager:

    #         with mock.patch("HPC_Drug.PDB.repair_pdb.pdbfixer.pdbfixer.PDBFixer", autospec = True) as mocked_PDBFixer:

    #             with mock.patch("HPC_Drug.PDB.repair_pdb.pdbfixer.simtk.openmm.app.pdbxfile.PDBxFile.writeFile", autospec = True) as mocked_pdbwriter:

    #                 mocked_context_manager.return_value.__enter__.return_value = "StringIO"

    #                 mocked_fixer = mock.Mock()

    #                 mocked_PDBFixer.return_value = mocked_fixer
            
    #                 #The function to test
    #                 output = pdbfixer.repair(input_file_name = "test_input.cif",
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
