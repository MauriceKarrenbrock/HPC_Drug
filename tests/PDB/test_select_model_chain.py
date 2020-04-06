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

from HPC_Drug.PDB import select_model_chain

class test_select_model_chain(unittest.TestCase):

    def test_with_correct_input(self):

        #I did a very bad job here, will do better sooner or later
        
        Protein = mock.MagicMock()

        Protein.protein_id = "id"

        Protein.pdb_file = "test_file"

        Protein.model = 0

        Protein.chain = "A"

        Protein.structure.return_value = { 0 : { "A" : Protein } }

        select_model_chain.select_model_chain(Protein = Protein)

        Protein.update_structure.assert_called_once()

        Protein.write.assert_called_once()