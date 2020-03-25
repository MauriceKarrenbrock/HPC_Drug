"""
Copyright (c) 2020-2020 Maurice Karrenbrock

This software is open-source and is distributed under the
GNU Affero General Public License v3 (agpl v3) license

A copy of the license must be included with any copy of the program or part of it
"""

import unittest
import Bio.PDB.Entity
from HPC_Drug import structure_informations
from HPC_Drug import structures


class test_StructureInformations(unittest.TestCase):

    def test_init_with_right_input(self):
        
        Protein = structures.Protein(protein_id="2f3z", filename = "tests/files4tests/2f3z.cif")

        test_class = structure_informations.StructureInformations(Protein = Protein)

        self.assertIsInstance(test_class.Protein, structures.Protein)

    def test_init_with_wrong_input(self):
        
        Protein = ("dummy", 1, 1.1)

        with self.assertRaises(TypeError):
            for i in Protein:
                test_class = structure_informations.StructureInformations(Protein = i)

    def test_get_structure_with_right_input(self):

        Protein = structures.Protein(protein_id="2f3z", filename = "tests/files4tests/2f3z.cif")

        test_class = structure_informations.StructureInformations(Protein = Protein)

        filename = (("tests/files4tests/2f3z.cif", "cif"), ("tests/files4tests/2f3z.pdb", "pdb"))

        for i in filename:
            self.assertIsInstance(test_class.get_structure(
                                                        filename = i[0], file_type = i[1]),
                                                        Bio.PDB.Entity.Entity)

    def test_get_structure_with_no_input_should_use_self_Protein(self):

        Protein = structures.Protein(protein_id="2f3z", filename = "tests/files4tests/2f3z.cif")

        test_class = structure_informations.StructureInformations(Protein = Protein)

        self.assertIsInstance(test_class.get_structure(),
                                                    Bio.PDB.Entity.Entity)

    def test_get_structure_with_wrong_filetype_type(self):

        Protein = structures.Protein(protein_id="2f3z", filename = "tests/files4tests/2f3z.cif")

        test_class = structure_informations.StructureInformations(Protein = Protein)

        filename = (("tests/files4tests/2f3z.cif", "WRONG"), ("tests/files4tests/2f3z.pdb", "WRONG"))

        with self.assertRaises(TypeError):
            for i in filename:
                self.assertIsInstance(test_class.get_structure(
                                                            filename = i[0], file_type = i[1]),
                                                            Bio.PDB.Entity.Entity)