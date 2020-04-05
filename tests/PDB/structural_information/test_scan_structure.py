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

from HPC_Drug.PDB.structural_information import scan_structure


class test_get_metal_binding_residues_with_no_header(unittest.TestCase):

    def test_with_wrong_structure(self):

        with self.assertRaises(TypeError):
        
            scan_structure.get_metal_binding_residues_with_no_header(structure = "WRONG",
                                                    protein_chain = 'A',
                                                    protein_model = 0)

        with self.assertRaises(AttributeError):
        
            scan_structure.get_metal_binding_residues_with_no_header(structure = "WRONG",
                                                    protein_chain = None,
                                                    protein_model = None)

    
    def test_with_no_model_chain_selection(self):

        with mock.patch("HPC_Drug.orient.Orient", autospec = True) as mocked_class:

            mocked_class.return_value.center_mass_distance.side_effect = (("DUMMY", "DUMMY", 0.), ("DUMMY", "DUMMY", 200.))

            #instantiating the residues
            residue_near_residue = mock.MagicMock()
            residue_far_residue = mock.MagicMock()
            residue_metal_residue = mock.MagicMock()

            #giving them resnames
            residue_near_residue.resname = "HIS"
            residue_far_residue.resname = "CYS"
            residue_metal_residue.resname = "ZN"

            #giving them resnums
            residue_near_residue.id = ("", 1, "")
            residue_far_residue.id = ("", 2, "")
            residue_metal_residue.id = ("h", 3, "")

            #instantiating some atoms
            near_atom = mock.MagicMock()
            far_atom = mock.MagicMock()
            zn_atom = mock.MagicMock()

            #giving the atoms names
            near_atom.name = "E2"
            far_atom.name = "FAR_AWAY"
            zn_atom.name = "ZN"

            #giving the atoms some coordinates
            near_atom.coord = (1., 0., 0.)
            far_atom.coord = (1.1, 0., 0.)
            zn_atom.coord = (0., 0., 0.)

            residue_near_residue.__iter__.return_value = [near_atom, far_atom]
            residue_metal_residue.__iter__.return_value = [zn_atom]

            residues = [residue_near_residue, residue_far_residue, residue_metal_residue]

            structure = mock.MagicMock()
            structure.return_value = structure
            structure.__iter__.return_value = residues
            structure.return_value.get_residues.return_value = residues

            expected_output = {1 : ["HIS", "E2", "ZN"]}

            output = scan_structure.get_metal_binding_residues_with_no_header(structure = structure,
                                                                            cutoff = 3.0,
                                                                            protein_chain = None,
                                                                            protein_model = None,
                                                                            COM_distance = 10.0,
                                                                            metals = ["ZN"])

            self.assertEqual(output, expected_output)


    def test_with_model_chain_selection(self):

        with mock.patch("HPC_Drug.orient.Orient", autospec = True) as mocked_class:

            mocked_class.return_value.center_mass_distance.side_effect = (("DUMMY", "DUMMY", 0.), ("DUMMY", "DUMMY", 200.))

            #instantiating the residues
            residue_near_residue = mock.MagicMock()
            residue_far_residue = mock.MagicMock()
            residue_metal_residue = mock.MagicMock()

            #giving them resnames
            residue_near_residue.resname = "HIS"
            residue_far_residue.resname = "CYS"
            residue_metal_residue.resname = "ZN"

            #giving them resnums
            residue_near_residue.id = ("", 1, "")
            residue_far_residue.id = ("", 2, "")
            residue_metal_residue.id = ("h", 3, "")

            #instantiating some atoms
            near_atom = mock.MagicMock()
            far_atom = mock.MagicMock()
            zn_atom = mock.MagicMock()

            #giving the atoms names
            near_atom.name = "E2"
            far_atom.name = "FAR_AWAY"
            zn_atom.name = "ZN"

            #giving the atoms some coordinates
            near_atom.coord = (1., 0., 0.)
            far_atom.coord = (1.1, 0., 0.)
            zn_atom.coord = (0., 0., 0.)

            residue_near_residue.__iter__.return_value = [near_atom, far_atom]
            residue_metal_residue.__iter__.return_value = [zn_atom]

            residues = [residue_near_residue, residue_far_residue, residue_metal_residue]

            structure = mock.MagicMock()
            structure.return_value = structure
            structure.__iter__.return_value = residues
            structure.return_value.get_residues.return_value = residues

            input_structure = { 0 : { 'A' : structure } }

            expected_output = {1 : ["HIS", "E2", "ZN"]}

            output = scan_structure.get_metal_binding_residues_with_no_header(structure = input_structure,
                                                                            cutoff = 3.0,
                                                                            protein_chain = 'A',
                                                                            protein_model = 0,
                                                                            COM_distance = 10.0,
                                                                            metals = ["ZN"])

            self.assertEqual(output, expected_output)



class test_get_disulf_bonds_with_no_header(unittest.TestCase):
    
    
    def test_with_wrong_structure(self):

        with self.assertRaises(TypeError):
        
            scan_structure.get_disulf_bonds_with_no_header(structure = "WRONG",
                                                    protein_chain = 'A',
                                                    protein_model = 0)

        with self.assertRaises(AttributeError):
        
            scan_structure.get_disulf_bonds_with_no_header(structure = "WRONG",
                                                    protein_chain = None,
                                                    protein_model = None)

    
    def test_with_no_model_chain_selection(self):
    
        #instantiating the residues
        residue_binding_1 = mock.MagicMock()
        residue_far_residue = mock.MagicMock()
        residue_binding_2 = mock.MagicMock()
        residue_not_CYS = mock.MagicMock()

        #giving them resnames
        residue_binding_1.resname = "CYS"
        residue_far_residue.resname = "CYS"
        residue_binding_2.resname = "CYS"
        residue_not_CYS.resname = "NOT_CYS"

        #giving them resnums
        residue_binding_1.id = ("", 1, "")
        residue_far_residue.id = ("", 2, "")
        residue_binding_2.id = ("", 3, "")
        residue_not_CYS.id = ("", 4, "")

        #instantiating some atoms
        SG_atom_binding_1 = mock.MagicMock()
        SG_atom_far_residue = mock.MagicMock()
        SG_atom_binding_2 = mock.MagicMock()

        #giving the atoms names
        SG_atom_binding_1.name = "SG"
        SG_atom_far_residue.name = "SG"
        SG_atom_binding_2.name = "SG"

        #giving the atoms some coordinates
        SG_atom_binding_1.coord = (1., 0., 0.)
        SG_atom_far_residue.coord = (100., 0., 0.)
        SG_atom_binding_2.coord = (0., 0., 0.)

        #the method to get the SG atom from the residue
        residue_binding_1.__getitem__.side_effect = { 'SG' : SG_atom_binding_1 }.__getitem__
        residue_far_residue.__getitem__.side_effect = { 'SG' : SG_atom_far_residue }.__getitem__
        residue_binding_2.__getitem__.side_effect = { 'SG' : SG_atom_binding_2 }.__getitem__

        residues = [residue_binding_1, residue_far_residue, residue_binding_2, residue_not_CYS]

        structure = mock.MagicMock()
        structure.return_value = structure
        structure.__iter__.return_value = residues
        structure.return_value.get_residues.return_value = residues

        expected_output = ( { 1 : ["CYS", "SG", "disulf"], 3 : ["CYS", "SG", "disulf"] } , [ ( 1 , 3 ) ] )

        output = scan_structure.get_disulf_bonds_with_no_header(structure = structure,
                                                                cutoff = 3.0,
                                                                protein_chain = None,
                                                                protein_model = None)

        self.assertEqual(output, expected_output)


    def test_with_model_chain_selection(self):

        #instantiating the residues
        residue_binding_1 = mock.MagicMock()
        residue_far_residue = mock.MagicMock()
        residue_binding_2 = mock.MagicMock()
        residue_not_CYS = mock.MagicMock()

        #giving them resnames
        residue_binding_1.resname = "CYS"
        residue_far_residue.resname = "CYS"
        residue_binding_2.resname = "CYS"
        residue_not_CYS.resname = "NOT_CYS"

        #giving them resnums
        residue_binding_1.id = ("", 1, "")
        residue_far_residue.id = ("", 2, "")
        residue_binding_2.id = ("", 3, "")
        residue_not_CYS.id = ("", 4, "")

        #instantiating some atoms
        SG_atom_binding_1 = mock.MagicMock()
        SG_atom_far_residue = mock.MagicMock()
        SG_atom_binding_2 = mock.MagicMock()

        #giving the atoms names
        SG_atom_binding_1.name = "SG"
        SG_atom_far_residue.name = "SG"
        SG_atom_binding_2.name = "SG"

        #giving the atoms some coordinates
        SG_atom_binding_1.coord = (1., 0., 0.)
        SG_atom_far_residue.coord = (100., 0., 0.)
        SG_atom_binding_2.coord = (0., 0., 0.)

        #the method to get the SG atom from the residue
        residue_binding_1.__getitem__.side_effect = { 'SG' : SG_atom_binding_1 }.__getitem__
        residue_far_residue.__getitem__.side_effect = { 'SG' : SG_atom_far_residue }.__getitem__
        residue_binding_2.__getitem__.side_effect = { 'SG' : SG_atom_binding_2 }.__getitem__

        residues = [residue_binding_1, residue_far_residue, residue_binding_2, residue_not_CYS]

        structure = mock.MagicMock()
        structure.return_value = structure
        structure.__iter__.return_value = residues
        structure.return_value.get_residues.return_value = residues

        expected_output = ( { 1 : ["CYS", "SG", "disulf"], 3 : ["CYS", "SG", "disulf"] } , [ ( 1 , 3 ) ] )

        input_structure = { 0 : { 'A' : structure } }

        output = scan_structure.get_disulf_bonds_with_no_header(structure = input_structure,
                                                                cutoff = 3.0,
                                                                protein_chain = 'A',
                                                                protein_model = 0)

        self.assertEqual(output, expected_output)



class test_get_organic_ligands_with_no_header(unittest.TestCase):
    
    def test_with_no_ligands(self):

        residue_1 = mock.Mock()
        residue_2 = mock.Mock()

        residue_1.resname = "ZN"
        residue_2.resname = "ZN"

        residue_1.id = ("h", 1, "")
        residue_2.id = ("h", 2, "")

        residues = [residue_1, residue_2]

        structure = mock.MagicMock()
        structure.return_value = structure
        structure.__iter__.return_value = residues

        output = scan_structure.get_organic_ligands_with_no_header(structure = structure,
                                                    protein_chain = None,
                                                    protein_model = None)

        self.assertEqual(output, None)

    def test_with_wrong_structure(self):

        with self.assertRaises(TypeError):
        
            scan_structure.get_organic_ligands_with_no_header(structure = "WRONG",
                                                    protein_chain = 'A',
                                                    protein_model = 0)

        with self.assertRaises(AttributeError):
        
            scan_structure.get_organic_ligands_with_no_header(structure = "WRONG",
                                                    protein_chain = None,
                                                    protein_model = None)

    def test_wit_no_chain_model_selection(self):

        residue_1 = mock.Mock()
        residue_2 = mock.Mock()

        residue_1.resname = "LIGAND"
        residue_2.resname = "ZN"

        residue_1.id = ("h", 1, "")
        residue_2.id = ("h", 2, "")

        residues = [residue_1, residue_2]

        structure = mock.MagicMock()
        structure.return_value = structure
        structure.__iter__.return_value = residues

        output = scan_structure.get_organic_ligands_with_no_header(structure = structure,
                                                    protein_chain = None,
                                                    protein_model = None,
                                                    metals = ["ZN"])

        self.assertEqual(output, [["LIGAND", 1]])


    def test_with_model_chain_selection(self):

        residue_1 = mock.Mock()
        residue_2 = mock.Mock()

        residue_1.resname = "LIGAND"
        residue_2.resname = "ZN"

        residue_1.id = ("h", 1, "")
        residue_2.id = ("h", 2, "")

        residues = [residue_1, residue_2]

        structure = mock.MagicMock()
        structure.return_value = structure
        structure.__iter__.return_value = residues

        input_structure = { 0 : { 'A' : structure } }

        output = scan_structure.get_organic_ligands_with_no_header(structure = input_structure,
                                                                protein_chain = 'A',
                                                                protein_model = 0,
                                                                metals = "ZN")

        self.assertEqual(output, [["LIGAND", 1]])



class test_get_metalbinding_disulf_ligands(unittest.TestCase):
    
    def test_with_correct_input(self):

        with mock.patch("HPC_Drug.PDB.structural_information.scan_structure.get_metal_binding_residues_with_no_header", return_value = {"metalc" : "metalc"}):

            with mock.patch("HPC_Drug.PDB.structural_information.scan_structure.get_disulf_bonds_with_no_header", return_value = ( {"disulf" : "disulf"}, [(1,2)] ) ):

                with mock.patch("HPC_Drug.PDB.structural_information.scan_structure.get_organic_ligands_with_no_header", return_value = [["LIGAND", 1]] ):

                    Protein = mock.Mock()

                    Protein, organic_ligands = scan_structure.get_metalbinding_disulf_ligands(Protein = Protein)

                    Protein.update_structure.assert_called_once_with(struct_type = "biopython")

                    self.assertEqual(Protein.substitutions_dict, { "metalc" : "metalc" , "disulf" : "disulf" })
                    self.assertEqual(Protein.sulf_bonds, [(1,2)])
                    self.assertEqual(organic_ligands, [["LIGAND", 1]])


                    