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

from HPC_Drug.structures import protein

class test_Protein(unittest.TestCase):

    def test_init(self):

        test_class = protein.Protein(protein_id = "protein_id",
                                pdb_file = None,
                                structure = "structure",
                                substitutions_dict = "substitutions_dict",
                                sulf_bonds = "sulf_bonds",
                                seqres = "seqres",
                                file_type = 'cif',
                                model = '0',
                                chain = 'a',
                                cys_dict = "cys_dict",
                                gro_file = "gro_file",
                                top_file = "top_file",
                                tpg_file = "tpg_file",
                                prm_file = "prm_file")

        self.assertEqual(test_class.protein_id, "protein_id")
        self.assertEqual(test_class.pdb_file, f"{test_class.protein_id}.{test_class.file_type}")
        self.assertEqual(test_class.structure, "structure")
        self.assertEqual(test_class.substitutions_dict, "substitutions_dict")
        self.assertEqual(test_class.sulf_bonds, "sulf_bonds")
        self.assertEqual(test_class.seqres, "seqres")
        self.assertEqual(test_class.file_type, "cif")
        self.assertEqual(test_class.model, 0)
        self.assertEqual(test_class.chain, 'A')
        self.assertEqual(test_class.cys_dict, "cys_dict")
        self.assertEqual(test_class.gro_file, "gro_file")
        self.assertEqual(test_class.top_file, "top_file")
        self.assertEqual(test_class.tpg_file, "tpg_file")
        self.assertEqual(test_class.prm_file, "prm_file")
        self.assertEqual(test_class._ligands, [])

    def test_init_raise(self):

        with self.assertRaises(ValueError):

            protein.Protein()

    def test_add_ligand(self):

        test_class = protein.Protein(protein_id = "test")

        test_class.add_ligand("test")

        self.assertEqual(test_class._ligands, ["test"])

    def test_get_ligand_list(self):

        test_class = protein.Protein(protein_id = "test")

        test_class._ligands = ["a", "b", "c"]

        output = test_class.get_ligand_list()

        self.assertIs(output, test_class._ligands)