######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

import pytest

import HPC_Drug.MD.openbabel_utils as _ob

class Test_convert_and_protonate_file_to_sdf():

    @pytest.mark.parametrize(
        'input_file',
        [('lig_no_H.pdb'),
         ('lig_no_H.sdf'),
         ('lig_no_H.mol2'),
         ('lig_no_H.smi')
         ])
    def test_with_ph7(self, get_data_dir, input_file, tmp_path):

        output_file = tmp_path / 'output.sdf'

        input_file = get_data_dir / input_file
        expected = get_data_dir / 'lig_yes_H.sdf'

        _ob.convert_and_protonate_file_to_sdf(input_file, output_file, ph=7, ligand_resname='UFV')

        with output_file.open() as out:
            with expected.open() as exp:
                assert out.readlines()[0] == exp.readlines()[0]
                assert out.readlines()[2:] == exp.readlines()[2:]


