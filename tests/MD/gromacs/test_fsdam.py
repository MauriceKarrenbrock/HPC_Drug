######################################################################################
# Copyright (c) 2020-2023 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

import mdtraj
import HPC_Drug.MD.gromacs.fsdam as _fsdam

class Test_create_complex_lig_creation_frames():

    def test_works(self, get_data_dir, tmp_path):
        num_files = _fsdam.create_complex_lig_creation_frames(apo_files=[str(get_data_dir / "apo.pdb")] * 3,
                                    holo_files=[str(get_data_dir / "holo.pdb")] * 3,
                                    ligand_in_vacuum_files=[str(get_data_dir / "lig_4_apo.pdb")] * 3,
                                    apo_files_top=str(get_data_dir / "apo.pdb"), 
                                    holo_files_top=str(get_data_dir / "holo.pdb"),
                                    ligand_in_vacuum_files_top=str(get_data_dir / "lig_4_apo.pdb"),
                                    lig_selection_string="resname LIG",
                                    output_dir=tmp_path,
                                    verbose=False)

        assert num_files == 3

        ref = mdtraj.load(str(get_data_dir / "apo_with_lig.pdb"))

        for gro in tmp_path.glob("*.pdb"):
            assert ref == mdtraj.load(str(gro))
                  