"""
Copyright (c) 2020-2020 Maurice Karrenbrock

This software is open-source and is distributed under the
GNU Affero General Public License v3 (agpl v3) license

A copy of the license must be included with any copy of the program or part of it
"""

import unittest
import os

from HPC_Drug.files_IO import write_on_files

class test_write_file(unittest.TestCase):

    def test_wrong_input(self):

        wrong_inputs = (1, 1.5)

        with self.assertRaises(TypeError):

            for i in wrong_inputs:
                
                write_on_files.write_file(lines = i, file_name = "dummy")

    
    def test_with_string(self):

        small_string = "string"
        string = f"{small_string}\n{small_string}"
        output_file = "tests/files4tests/test_with_string.txt"

        try:
            write_on_files.write_file(lines = string, file_name = output_file)

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
                write_on_files.write_file(lines = input_lists[i], file_name = output_file)

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
                write_on_files.write_file(lines = input_tuples[i], file_name = output_file)

                with open(output_file, 'r') as f:
                    lines = f.readlines()

                self.assertEqual(lines, expected_ouputs[i])
        finally:

            os.remove(output_file)


def _create_input_file():
    """private"""

    file_name = 'tests/files4tests/file_to_append.txt'
    string = "aaa"

    write_on_files.write_file(lines = string, file_name = file_name)

    return file_name, string

class test_append_file(unittest.TestCase):

    def test_wrong_input(self):

        wrong_inputs = (1, 1.5)

        with self.assertRaises(TypeError):

            for i in wrong_inputs:
                
                write_on_files.append_file(lines = i, file_name = "dummy")

    
    def test_with_string(self):

        input_file, existing_string = _create_input_file()

        small_string = "string"
        string = f"{small_string}\n{small_string}"

        try:
            write_on_files.append_file(lines = string, file_name = input_file)

            with open(input_file, 'r') as f:
                lines = f.readlines()

            self.assertEqual(lines, [existing_string + small_string + '\n', small_string])
        finally:
            
            os.remove(input_file)

    def test_with_list(self):

        input_lists = (["string", "string"], ["\nstring\n", "string\n"])
        

        for i in range(len(input_lists)):

            input_file, existing_string = _create_input_file()
            expected_ouputs = ([existing_string + "stringstring"], [existing_string + '\n', "string\n", "string\n"])

            try:

                write_on_files.append_file(lines = input_lists[i], file_name = input_file)

                with open(input_file, 'r') as f:
                    lines = f.readlines()

                self.assertEqual(lines, expected_ouputs[i])

            finally:

                os.remove(input_file)

    def test_with_tuple(self):

        input_tuples = (("string", "string"), ("\nstring\n", "string\n"))
        

        for i in range(len(input_tuples)):

            input_file, existing_string = _create_input_file()
            expected_ouputs = ([existing_string + "stringstring"], [existing_string + '\n', "string\n", "string\n"])

            try:

            
                write_on_files.append_file(lines = input_tuples[i], file_name = input_file)

                with open(input_file, 'r') as f:
                    lines = f.readlines()

                self.assertEqual(lines, expected_ouputs[i])

            finally:

                os.remove(input_file)

