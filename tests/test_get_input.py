import unittest
import unittest.mock
import importlib_resources

import HPC_Drug.get_input


class test_GetInput(unittest.TestCase):

    def test_init_function_with_right_inputs(self):

        input_class = HPC_Drug.get_input.GetInput(1, 'a', 1.0, k = 1, r = 2)

        self.assertEqual(input_class.input, (1, 'a', 1.0))
        self.assertEqual(input_class.dictionary, {'k':1, 'r':2})


class test_GetFile(unittest.TestCase):

    def test_init_function_with_right_filename(self):
        filename = 'test_filename'
        file_class = HPC_Drug.get_input.GetFile(filename)
        self.assertEqual(file_class.filename, filename)

    def test_init_function_with_wrong_filetype(self):
        filename = (3.4, 3, None)
        for name in filename:
            with self.assertRaises(TypeError):
                file_class = HPC_Drug.get_input.GetFile(filename = name)

    def test_init_function_with_no_file(self):
        with self.assertRaises(TypeError):
            file_class = HPC_Drug.get_input.GetFile()


class ParseInputFromFile(unittest.TestCase):

    def test_create_input_dict(self):
        
        #dummy filename to make the __init__ function happy
        filename = "tests/input_correct_4tests.txt"

        test_class = HPC_Drug.get_input.ParseInputFromFile(filename)
        dictionary = test_class._create_input_dict()

        self.assertEqual(set(dictionary.keys()), set(test_class.possible_keys))
        for value in dictionary.values():
            self.assertEqual(value, None)

    
    def test_refine_input_correct_input(self):

        filename = "tests/input_correct_4tests.txt"

        input_dict = {'Protein_model' : '0',
                    'Protein_chain' : None,
                    'ph' : '7.0',
                    'ligand_elaboration_program' : None,
                    'ligand_elaboration_program_path' : None,
                    'protein_tpg_file' : None,
                    'protein_prm_file' : None,
                    'solvent_pdb' : None}

        with importlib_resources.path('HPC_Drug.lib', 'amber99sb-ildn.tpg') as tpg:
            with importlib_resources.path('HPC_Drug.lib', 'amber99sb-ildn.prm') as prm:
                with importlib_resources.path('HPC_Drug.lib', 'water.pdb') as solv:
                    expected_dict = {'Protein_model' : 0,
                                    'Protein_chain' : 'A',
                                    'ph' : 7.0,
                                    'ligand_elaboration_program' : 'primadorac',
                                    'ligand_elaboration_program_path' :  '~/ORAC/trunk/tools/primadorac/primadorac.bash',
                                    'protein_tpg_file' : tpg,
                                    'protein_prm_file' : prm,
                                    'solvent_pdb' : solv}
        
        test_class = HPC_Drug.get_input.ParseInputFromFile(filename)

        output_dict = test_class._refine_input(input_dict)

        self.assertEqual(expected_dict, output_dict)

    def test_refine_input_wrong_input_type(self):

        filename = "tests/input_correct_4tests.txt"

        wrong_input = (77.3, 'a', None)

        test_class = HPC_Drug.get_input.ParseInputFromFile(filename)

        for wrong in wrong_input:
            with self.assertRaises(TypeError):
                test_class._refine_input(wrong)

    def test_refine_input_missing_dict_keys(self):

        filename = "tests/input_correct_4tests.txt"

        wrong_input_dict = {'Protein_chain' : None,
                            'ph' : '7.0',
                            'ligand_elaboration_program_path' : None,
                            'protein_tpg_file' : None,
                            'protein_prm_file' : None,
                            'solvent_pdb' : None}

        test_class = HPC_Drug.get_input.ParseInputFromFile(filename)

        with self.assertRaises(KeyError):
            test_class._refine_input(wrong_input_dict)


    def test_read_input_with_correct_input(self):
        """Shall be done better, I didn't mock the private methods"""

        filename = "tests/input_correct_4tests.txt"
        test_class = HPC_Drug.get_input.ParseInputFromFile(filename)

        dictionary = test_class.read_input()

        self.assertEqual(set(dictionary.keys()), set(test_class.possible_keys))

        #da rifare daccapo
        self.assertEqual(dictionary, {'protein' : '1df8',
                                    'ligand' : None,
                                    'protein_filetype' : 'cif',
                                    'ligand_elaboration_program' : 'primadorac',
                                    'ligand_elaboration_program_path' : '~/ORAC/trunk/tools/primadorac/primadorac.bash',
                                    'local' : 'no',
                                    'filepath' : None,
                                    'Protein_model' : 0,
                                    'Protein_chain' : 'A',
                                    'ligand_in_protein' : 'yes',
                                    'ph' : 7.0,
                                    'repairing_method' : 'pdbfixer',
                                    'MD_program' : 'orac',
                                    'MD_program_path' : '~/ORAC/trunk/src/GNU-FFTW-OMP/orac',
                                    'protein_prm_file' : 'amber99sb-ildn.prm',
                                    'protein_tpg_file' : 'amber99sb-ildn.tpg',
                                    'solvent_pdb' : 'water.pdb'})


    def test_read_input_with_wrong_input_key(self):
        """Shall be done better, I didn't mock the private methods"""
        
        filename = "tests/input_wrong.txt"

        with self.assertRaises(ValueError):
            test_class = HPC_Drug.get_input.ParseInputFromFile(filename)   

if __name__ == '__main__':
    unittest.main()