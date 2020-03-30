"""
Copyright (c) 2020-2020 Maurice Karrenbrock

This software is open-source and is distributed under the
GNU Affero General Public License v3 (agpl v3) license

A copy of the license must be included with any copy of the program or part of it
"""

import unittest
import unittest.mock

import HPC_Drug.orient
import HPC_Drug.structures.ligand
import HPC_Drug.structures.protein


class TestOrient(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        Protein = HPC_Drug.structures.protein.Protein(protein_id = '2gz7',
                pdb_file = 'tests/files4tests/2gz7.cif',
                structure = None,
                substitutions_dict = None,
                sulf_bonds = None,
                seqres = None,
                file_type = 'cif',
                model = 0,
                chain = 'A',
                gro_file = None,
                top_file = None)

        Ligand = HPC_Drug.structures.ligand.Ligand(resname = 'D3F',
                pdb_file = None,
                structure = None,
                file_type = 'pdb',
                tpg_file = None,
                prm_file = None,
                resnum = 307,
                itp_file = None)

        cls.orient_class = HPC_Drug.orient.Orient(Protein = Protein,
                                                Ligand = Ligand)

    def test_get_hot_residues_for_rem_with_cif(self):

        right_list = [['HIS', 41], ['HIS', 164], ['GLN', 192]]
        residues_list = self.orient_class.get_hot_residues_for_rem(cutoff = 3.0, residue_dist = 8.0)

        #sort them in order not to worry about order
        self.assertEqual(right_list.sort(), residues_list.sort())

