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

from HPC_Drug.PDB.organic_ligand import get_ligand_topology

class test_get_ligand_topology(unittest.TestCase):

    def test_raise_NotImplementedError(self):

        with self.assertRaises(NotImplementedError):

            with mock.patch("HPC_Drug.auxiliary_functions.path.absolute_programpath"):

                Protein = mock.Mock()

                Protein.get_ligand_list.return_value = ["dummy_ligand"]

                get_ligand_topology.get_topology(Protein = Protein, program_path = "test_path", tool = "WRONG_TOOL", ph = 7.0)

    
    def test_empty_ligand(self):

        Protein = mock.Mock()

        Protein.return_value = Protein

        Protein.get_ligand_list.return_value = []

        with mock.patch("HPC_Drug.PDB.organic_ligand.primadorac.Primadorac") as mocked_primadorac:

            with mock.patch("HPC_Drug.auxiliary_functions.path.absolute_programpath") as mocked_abs_programpath:

                output = get_ligand_topology.get_topology(Protein = Protein, program_path = "test_path", tool = "primadorac", ph = 7.0)

                mocked_primadorac.assert_not_called()

                mocked_abs_programpath.assert_not_called()

                self.assertEqual(output.get_ligand_list(), [])

    def test_None_ligand(self):

        Protein = mock.Mock()

        Protein.return_value = Protein

        Protein.get_ligand_list.return_value = None

        with mock.patch("HPC_Drug.PDB.organic_ligand.primadorac.Primadorac") as mocked_primadorac:

            with mock.patch("HPC_Drug.auxiliary_functions.path.absolute_programpath") as mocked_abs_programpath:

                output = get_ligand_topology.get_topology(Protein = Protein, program_path = "test_path", tool = "primadorac", ph = 7.0)

                mocked_primadorac.assert_not_called()

                mocked_abs_programpath.assert_not_called()

                self.assertEqual(output.get_ligand_list(), None)

    
    def test_tool_primadorac(self):

        Protein = mock.Mock()

        Protein.return_value = Protein

        Protein.get_ligand_list.return_value = ["test_ligand"]

        with mock.patch("HPC_Drug.PDB.organic_ligand.primadorac.Primadorac") as mocked_primadorac:

            with mock.patch("HPC_Drug.auxiliary_functions.path.absolute_programpath", return_value = "test_path") as mocked_abs_programpath:

                get_ligand_topology.get_topology(Protein = Protein, program_path = "test_path", tool = "primadorac", ph = 7.0)

                mocked_primadorac.assert_called_once_with(Protein = Protein, primadorac_path = "test_path", ph = 7.0)

                mocked_primadorac.return_value.execute.assert_called_once()

                mocked_abs_programpath.assert_called_once()


