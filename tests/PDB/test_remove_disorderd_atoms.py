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

from HPC_Drug.PDB import remove_disordered_atoms

class test_remove_disordered_atoms(unittest.TestCase):

    def test_cif(self):

        with mock.patch("HPC_Drug.PDB.remove_disordered_atoms.Bio.PDB.MMCIFParser", autospec = True) as mocked_parser:

            with mock.patch("HPC_Drug.PDB.remove_disordered_atoms.Bio.PDB.MMCIFIO") as mocked_IO:

                Protein = mock.Mock()

                Protein.file_type = "cif"

                Protein.pdb_file = "test_file"

                Protein.structure = "test_structure"

                Protein.protein_id = "id"

                remove_disordered_atoms.remove_disordered_atoms(Protein = Protein)

                mocked_parser.assert_called_once()
                mocked_parser.return_value.assert_has_calls([mock.call.get_structure("id", "test_file")])

                mocked_IO.assert_called_once()
                mocked_IO.return_value.set_structure.assert_called_once()
                mocked_IO.return_value.save.assert_called_once()


    def test_pdb(self):

        with mock.patch("HPC_Drug.PDB.remove_disordered_atoms.Bio.PDB.PDBParser", autospec = True) as mocked_parser:

            with mock.patch("HPC_Drug.PDB.remove_disordered_atoms.Bio.PDB.PDBIO") as mocked_IO:

                Protein = mock.Mock()

                Protein.file_type = "pdb"

                Protein.pdb_file = "test_file"

                Protein.structure = "test_structure"

                Protein.protein_id = "id"

                remove_disordered_atoms.remove_disordered_atoms(Protein = Protein)

                mocked_parser.assert_called_once()
                mocked_parser.return_value.assert_has_calls([mock.call.get_structure("id", "test_file")])

                mocked_IO.assert_called_once()
                mocked_IO.return_value.set_structure.assert_called_once()
                mocked_IO.return_value.save.assert_called_once()

    def test_raise(self):

        with self.assertRaises(TypeError):

            Protein = mock.Mock()

            Protein.file_type = "WORNG_TYPE"

            remove_disordered_atoms.remove_disordered_atoms(Protein = Protein)