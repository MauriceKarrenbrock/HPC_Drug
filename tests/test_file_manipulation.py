######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################


import unittest

from HPC_Drug import file_manipulation
from HPC_Drug import pipeline_functions
from HPC_Drug.structures import ligand
from HPC_Drug.structures import protein
from HPC_Drug import important_lists

import Bio.PDB
import Bio.PDB.MMCIF2Dict
import os
import prody


class test_FileCruncer(unittest.TestCase):

    def test_init_function(self):
        with self.assertRaises(NotImplementedError):
            test_class = file_manipulation.FileCruncer()

class test_ProteinCruncer(unittest.TestCase):
    def test_ProteinCruncer(self):
        a = file_manipulation.ProteinCruncer('pdb')
        self.assertIsInstance(a, file_manipulation.PDBCruncer)

        a = file_manipulation.ProteinCruncer('cif')
        self.assertIsInstance(a, file_manipulation.MMCIFCruncer)

        with self.assertRaises(NotImplementedError):
            a = file_manipulation.ProteinCruncer('dummy')

class test_PDBCruncer(unittest.TestCase):
    
    def test_parse_must_return_atomgroup(self):

        pdb_file = "tests/files4tests/2f3z.pdb"
        protein_id = "2f3z"

        test_class = file_manipulation.PDBCruncer()

        self.assertIsInstance(test_class.parse(protein_id, pdb_file), prody.AtomGroup)
        self.assertIsInstance(test_class.parse(protein_id), prody.AtomGroup)

    #Too slow, because it connects to the pdb database
    # def test_parse_with_wrong_input_file(self):
    #     pdb_file = "tests/files4tests/dummy.pdb"
    #     protein_id = "dummy"

    #     test_class = file_manipulation.PDBCruncer()

    #     with self.assertRaises(OSError):
    #         test_class.parse(protein_id, pdb_file)

    #Too slow, because it connects to the pdb database
    # def test_parse_with_wrong_input_protein_id(self):
    #     protein_id = "dummy"

    #     test_class = file_manipulation.PDBCruncer()

    #     with self.assertRaises(OSError):
    #         test_class.parse(protein_id)

    def test_parse_with_wrong_input_type(self):
        protein_id = 4
        pdb_file = 2.2

        test_class = file_manipulation.PDBCruncer()

        with self.assertRaises(TypeError):
            test_class.parse(protein_id)
            test_class.parse(protein_id, pdb_file)

    def test_get_protein(self):
        #to do
        pass

    def test_get_ligand(self):
        #to do
        pass

    def test_get_ligand_with_no_ligands(self):
        
        test_class = file_manipulation.PDBCruncer()

        self.assertIsNone(test_class.get_ligand(
                                            "dummy_protein_id",
                                            "dummy_filename",
                                            None,
                                            prody.AtomGroup(),
                                            None))
    
    def test_get_ligand_with_wrong_structure_type(self):
        
        test_class = file_manipulation.PDBCruncer()

        with self.assertRaises(TypeError):
            test_class.get_ligand(
                                "dummy_protein_id",
                                "dummy_filename",
                                None,
                                "WRONG_STRUCTURE_TYPE",
                                None)

    def test_write_PDB(self):
        #to do
        pass
    




class test_MMCIFCruncer(unittest.TestCase):
    pass


class test_PDBRepair(unittest.TestCase):
    pass

class test_SubstitutionParser(unittest.TestCase):
    pass

class test_no_header_functions(unittest.TestCase):

    def test_get_metal_binding_residues_with_no_header(self):

        right_substitution_dict = {
            '179' : ['HIS', 'ND1', 'ZN'],
            '176' : ['CYS', 'SG', 'ZN'],
            '238' : ['CYS', 'SG', 'ZN'],
            '242' : ['CYS', 'SG', 'ZN']
        }

        output_dict = file_manipulation.get_metal_binding_residues_with_no_header(protein_id = '5aol',
                                                                                pdb_file = 'tests/files4tests/5aol.pdb',
                                                                                mmcif_file = None,
                                                                                cutoff = 3.0,
                                                                                substitutions_dict = {},
                                                                                protein_chain = 'A',
                                                                                protein_model = 0,
                                                                                COM_distance = 7.0)

        self.assertEqual(output_dict, right_substitution_dict)

    
    def test_get_disulf_bonds_with_no_header(self):

        right_substitution_dict = {
            '26' : ['CYS', 'SG', 'disulf'],
            '84' : ['CYS', 'SG', 'disulf'],
            '40' : ['CYS', 'SG', 'disulf'],
            '95' : ['CYS', 'SG', 'disulf'],
            '58' : ['CYS', 'SG', 'disulf'],
            '110' : ['CYS', 'SG', 'disulf'],
            '65' : ['CYS', 'SG', 'disulf'],
            '72' : ['CYS', 'SG', 'disulf']
        }
        right_sulf_bond = [('26', '84'), ('40', '95'), ('58', '110'), ('65', '72')]

        output_dict, output_sulf_bond = file_manipulation.get_disulf_bonds_with_no_header(protein_id = '3dxg',
                                                                                        pdb_file = 'tests/files4tests/3dxg.pdb',
                                                                                        mmcif_file = None,
                                                                                        cutoff = 3.0,
                                                                                        substitutions_dict = {},
                                                                                        sulf_bonds = [],
                                                                                        protein_chain = 'A',
                                                                                        protein_model = 0)

        self.assertEqual(output_dict, right_substitution_dict)
        #I sort the lists in order not to worry about order
        self.assertEqual(output_sulf_bond.sort(), right_sulf_bond.sort())

    
    def test_get_organic_ligands_with_no_header(self):


        right_ligand = [['U5P', '125'], ['U5P', '126']]


        output_ligand = file_manipulation.get_organic_ligands_with_no_header(protein_id = '3dxg',
                                                                            pdb_file = 'tests/files4tests/3dxg.pdb',
                                                                            mmcif_file = None,
                                                                            protein_chain = 'A',
                                                                            protein_model = 0)
        
        #I sort the lists in order not to worry about order
        self.assertEqual(output_ligand.sort(), right_ligand.sort())


        