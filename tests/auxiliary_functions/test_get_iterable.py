"""
Copyright (c) 2020-2020 Maurice Karrenbrock

This software is open-source and is distributed under the
GNU Affero General Public License v3 (agpl v3) license

A copy of the license must be included with any copy of the program or part of it
"""

import unittest

from HPC_Drug.auxiliary_functions import get_iterable

class test_get_iterable(unittest.TestCase):

    def test_with_string(self):

        output = get_iterable.get_iterable("test")

        self.assertEqual(output, ("test",) )

    def test_with_integer(self):

        output = get_iterable.get_iterable(1)

        self.assertEqual(output, (1,) )

    def test_with_list(self):

        output = get_iterable.get_iterable([1, "1"])

        self.assertEqual(output, [1, "1"] )