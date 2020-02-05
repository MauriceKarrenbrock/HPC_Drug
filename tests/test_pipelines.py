import unittest
import unittest.mock

import HPC_Drug.pipelines

#to do
class test_Pipeline(unittest.TestCase):

    def test_download_exception_wrong_local_input(self):

        wrong_names = ('ggg', 1, 1., ['a', 'yes', 'no'])

        for name in wrong_names:
            with self.assertRaises(ValueError):

                p = HPC_Drug.pipelines.Pipeline(local = name, filepath = 'dummy')
                file_class = p.download()