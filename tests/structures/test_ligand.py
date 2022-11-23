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

from HPC_Drug.structures import ligand

class test_Ligand(unittest.TestCase):

    def test_init(self):

        test_class = ligand.Ligand(resname = "resname",
                                file_type = 'pdb',
                                pdb_file = None,
                                structure = "structure",
                                resnum = '1',
                                itp_file = "itp_file",
                                gro_file = "gro_file",
                                top_file = "top_file",
                                tpg_file = "tpg_file",
                                prm_file = "prm_file")

        self.assertEqual(test_class.resname, "RESNAME")
        self.assertEqual(test_class.file_type, "pdb")
        self.assertEqual(test_class.pdb_file, f"{test_class.resname}_lgand.{test_class.file_type}")
        self.assertEqual(test_class.structure, "structure")
        self.assertEqual(test_class.resnum, 1)
        self.assertEqual(test_class.itp_file, "itp_file")
        self.assertEqual(test_class.gro_file, "gro_file")
        self.assertEqual(test_class.top_file, "top_file")
        self.assertEqual(test_class.tpg_file, "tpg_file")
        self.assertEqual(test_class.prm_file, "prm_file")
        
        
