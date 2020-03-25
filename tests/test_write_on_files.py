import unittest
import os

from HPC_Drug.files_IO import write_on_files

class test_write_file(unittest.TestCase):

    def test_wrong_input(self):

        wrong_inputs = (1, 1.5)

        with self.assertRaises(TypeError):

            for i in wrong_inputs:
                
                write_on_files.write_file(lines = i, file = "dummy")

    
    def test_with_string(self):

        small_string = "string"
        string = f"{small_string}\n{small_string}"
        output_file = "tests/files4tests/test_with_string.txt"

        try:
            write_on_files.write_file(lines = string, file = output_file)

            with open(output_file, 'r') as f:
                lines = f.readlines()

            self.assertEqual(lines, [small_string + '\n', small_string])
        finally:

            os.remove(output_file)

    def test_with_list(self):

        input_lists = (["string", "string"], ["string\n", "string\n"])
        expected_ouputs = (["stringstring"], ["string\n", "string\n"])

        output_file = "tests/files4tests/test_with_list.txt"

        try:

            for i in range(len(input_lists)):
                write_on_files.write_file(lines = input_lists[i], file = output_file)

                with open(output_file, 'r') as f:
                    lines = f.readlines()

                self.assertEqual(lines, expected_ouputs[i])
        finally:

            os.remove(output_file)

    def test_with_tuple(self):

        input_tuples = (("string", "string"), ("string\n", "string\n"))
        expected_ouputs = (["stringstring"], ["string\n", "string\n"])

        output_file = "tests/files4tests/test_with_tuple.txt"

        try:

            for i in range(len(input_tuples)):
                write_on_files.write_file(lines = input_tuples[i], file = output_file)

                with open(output_file, 'r') as f:
                    lines = f.readlines()

                self.assertEqual(lines, expected_ouputs[i])
        finally:

            os.remove(output_file)
