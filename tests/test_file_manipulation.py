import unittest

from HPC_Drug import file_manipulation
from HPC_Drug import pipeline_functions
from HPC_Drug import structures
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

        pdb_file = "tests/2f3z.pdb"
        protein_id = "2f3z"

        test_class = file_manipulation.PDBCruncer()

        self.assertIsInstance(test_class.parse(protein_id, pdb_file), prody.AtomGroup)
        self.assertIsInstance(test_class.parse(protein_id), prody.AtomGroup)

    #Too slow, because it connects to the pdb database
    # def test_parse_with_wrong_input_file(self):
    #     pdb_file = "tests/dummy.pdb"
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