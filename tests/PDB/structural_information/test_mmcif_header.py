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

from HPC_Drug.PDB.structural_information import mmcif_header


class test_get_metal_binding_residues(unittest.TestCase):

    def test_empty_dictionary_input(self):

        output = mmcif_header.get_metal_binding_residues(mmcif2dict = {})

        self.assertEqual(output, {})


    def test_not_dictionay_input(self):

        wrong_inputs = ("aaa", 1, 1., [1,1], (1,1))

        for i in wrong_inputs:

            with self.assertRaises(AttributeError):

                mmcif_header.get_metal_binding_residues(mmcif2dict = i)

    def test_with_valid_dictionary(self):

        input_dict = {
            '_struct_conn.conn_type_id' : ['metalc', 'metalc', 'DUMMY'],
            '_struct_conn.ptnr1_label_comp_id' : ['HIS', 'ZN', 'DUMMY'],
            '_struct_conn.ptnr2_label_comp_id' : ['ZN', 'HIS', 'DUMMY'],
            '_struct_conn.ptnr1_label_atom_id' : ['E2', 'ZN', 'DUMMY'],
            '_struct_conn.ptnr2_label_atom_id' : ['ZN', 'E2', 'DUMMY'],
            '_struct_conn.ptnr1_auth_seq_id' : ['1', '2', '3'],
            '_struct_conn.ptnr2_auth_seq_id' : ['4', '5', '6']
        }

        expected_output_dict = {
            1 : ('HIS', 'E2', 'ZN'),
            5 : ('HIS', 'E2', 'ZN')
        }

        output = mmcif_header.get_metal_binding_residues(mmcif2dict = input_dict, metals = ["ZN"])

        self.assertEqual(output, expected_output_dict)

    def test_not_implemented_metal(self):

        with mock.patch("HPC_Drug.PDB.structural_information.mmcif_header.print") as mocked_function:

            input_dict = {
                '_struct_conn.conn_type_id' : 'metalc',
                '_struct_conn.ptnr1_label_comp_id' : 'HIS',
                '_struct_conn.ptnr2_label_comp_id' : 'UNKNOWN',
                '_struct_conn.ptnr1_label_atom_id' : 'E2',
                '_struct_conn.ptnr2_label_atom_id' : 'UNKNOWN',
                '_struct_conn.ptnr1_auth_seq_id' : '1',
                '_struct_conn.ptnr2_auth_seq_id' : '4'
            }

            expected_output_dict = {}

            output = mmcif_header.get_metal_binding_residues(mmcif2dict = input_dict, metals = ["ZN"])

            mocked_function.assert_called_once()

            self.assertEqual(output, expected_output_dict)


class test_get_disulf_bonds(unittest.TestCase):

    def test_empty_dictionary_input(self):

        output = mmcif_header.get_disulf_bonds(mmcif2dict = {})

        self.assertEqual(output, ({}, []) )
 
    def test_with_valid_dictionary(self):

        input_dict = {
            '_struct_conn.conn_type_id' : ['disulf', 'disulf', 'DUMMY'],
            '_struct_conn.ptnr1_label_comp_id' : ['CYS', 'CYS', 'DUMMY'],
            '_struct_conn.ptnr2_label_comp_id' : ['CYS', 'CYS', 'DUMMY'],
            '_struct_conn.ptnr1_label_atom_id' : ['SG', 'SG', 'DUMMY'],
            '_struct_conn.ptnr2_label_atom_id' : ['SG', 'SG', 'DUMMY'],
            '_struct_conn.ptnr1_auth_seq_id' : ['1', '2', '3'],
            '_struct_conn.ptnr2_auth_seq_id' : ['2', '1', '4']
        }

        expected_output_dict = {
            1 : ('CYS', 'SG', 'disulf'),
            2 : ('CYS', 'SG', 'disulf')
        }

        expected_output_list = [(1, 2)]

        output_dict, output_list = mmcif_header.get_disulf_bonds(mmcif2dict = input_dict)

        self.assertEqual(output_dict, expected_output_dict)
        self.assertEqual(output_list, expected_output_list)


class test_get_organic_ligands(unittest.TestCase):

    def test_empty_dictionary_input(self):

        output = mmcif_header.get_organic_ligands(mmcif2dict = {})

        self.assertEqual(output, [])


    def test_not_dictionay_input(self):

        wrong_inputs = ("aaa", 1, 1., [1,1], (1,1))

        for i in wrong_inputs:

            with self.assertRaises(AttributeError):

                mmcif_header.get_organic_ligands(mmcif2dict = i)

    def test_no_chain_selection(self):

        input_dict = {
            '_struct_site.details' : ["LIGAND_A A DUMMY", "LIGAND_B B DUMMY", "LIGAND_A A DUMMY"]
        }

        expected_output = ["LIGAND_A", "LIGAND_B"]

        output = mmcif_header.get_organic_ligands(mmcif2dict = input_dict, protein_chain = None)

        self.assertEqual(sorted(output), sorted(expected_output))

    def test_with_chain_selection(self):

        input_dict = {
            '_struct_site.details' : ["LIGAND_A A DUMMY", "LIGAND_B B DUMMY", "LIGAND_A A DUMMY"]
        }

        expected_output = ["LIGAND_A"]

        output = mmcif_header.get_organic_ligands(mmcif2dict = input_dict, protein_chain = 'A')

        self.assertEqual(output, expected_output)


class test_get_ligand_resnum(unittest.TestCase):

    def test_with_no_ligands(self):

        self.assertEqual(mmcif_header.get_ligand_resnum(structure = "test", ligand_resnames = []), None)
        self.assertEqual(mmcif_header.get_ligand_resnum(structure = "test", ligand_resnames = None), None)

    def test_with_wrong_structure_type(self):

        with self.assertRaises(TypeError):

            mmcif_header.get_ligand_resnum(structure = "test_wrong", ligand_resnames = ["aa", "bb"])

    def test_with_no_chain_selection(self):

        input_ligand_resnames = ["LIGAND"]

        residue_1 = mock.Mock()
        residue_2 = mock.Mock()

        residue_1.resname = "LIGAND"
        residue_2.resname = "NOT_LIGAND"

        residue_1.id = ("", 1, "")
        residue_2.id = ("", 2, "")

        residues = [residue_1, residue_2]

        structure = mock.Mock()
        structure.return_value = structure
        structure.return_value.get_residues.return_value = residues

        output = mmcif_header.get_ligand_resnum(structure = structure,
                                            ligand_resnames = input_ligand_resnames,
                                            protein_chain = None,
                                            protein_model = None)

        self.assertEqual(output, [["LIGAND", 1]])

    def test_with_chain_model_selection(self):

        input_ligand_resnames = ["LIGAND"]

        residue_1 = mock.Mock()
        residue_2 = mock.Mock()

        residue_1.resname = "LIGAND"
        residue_2.resname = "NOT_LIGAND"

        residue_1.id = ("", 1, "")
        residue_2.id = ("", 2, "")

        residues = [residue_1, residue_2]

        structure = mock.Mock()
        structure.return_value = structure
        structure.return_value.get_residues.return_value = residues

        input_structure = {0 : {'A' : structure}}

        output = mmcif_header.get_ligand_resnum(structure = input_structure,
                                                ligand_resnames = input_ligand_resnames,
                                                protein_chain = 'A',
                                                protein_model = 0)

        self.assertEqual(output, [["LIGAND", 1]])

class test_get_metalbinding_disulf_ligands(unittest.TestCase):

    def test_with_correct_input(self):

        with mock.patch("HPC_Drug.PDB.structural_information.mmcif_header.get_metal_binding_residues", return_value = {"metalc" : "metalc"}):

            with mock.patch("HPC_Drug.PDB.structural_information.mmcif_header.get_disulf_bonds", return_value = ({'disulf':'disulf'}, [(1,2)])):

                with mock.patch("HPC_Drug.PDB.structural_information.mmcif_header.get_organic_ligands", return_value = ["ligand"]):

                    with mock.patch("HPC_Drug.PDB.structural_information.mmcif_header.get_ligand_resnum", return_value = [["ligand", 1]]):
                        
                        Protein = mock.Mock()
                        Protein.return_value.update_structure.return_value = "structure"

                        Protein, ligand_list = mmcif_header.get_metalbinding_disulf_ligands(Protein = Protein)

                        self.assertEqual(Protein.sulf_bonds, [(1,2)])
                        self.assertEqual(Protein.substitutions_dict, {"metalc" : "metalc", "disulf" : "disulf"})
                        self.assertEqual(ligand_list, [["ligand", 1]])