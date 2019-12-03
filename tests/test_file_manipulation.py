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
    pass

class test_MMCIFCruncer(unittest.TestCase):
    pass


class test_PDBRepair(unittest.TestCase):
    pass

class test_SubstitutionParser(unittest.TestCase):
    pass