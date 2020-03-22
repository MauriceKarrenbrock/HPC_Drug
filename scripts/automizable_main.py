"""
Copyright (c) 2020-2020 Maurice Karrenbrock

This software is open-source and is distributed under the
GNU Affero General Public License v3 (agpl v3) license

A copy of the license must be included with any copy of the program or part of it
"""

#This is an example of automated program that executes the given pipeline for any given protein id in a file
#that has a protein id per line
#es:
#1df8
#2gz7
#...

# python automizable_main.py input_file.txt

import HPC_Drug
import HPC_Drug.pipelines

import subprocess
import os
import sys
import contextlib

input_file_name = sys.argv[1]

with open(input_file_name, 'r') as f:
    lines = f.readlines()

for line in lines:
    protein_id = line.strip().lower()

    try:

        p = HPC_Drug.pipelines.NoLigand_Pipeline(protein = protein_id,
                                                protein_filetype = 'cif',
                                                local = 'no',
                                                filepath = None,
                                                ligand = None,
                                                ligand_elaboration_program = 'primadorac',
                                                ligand_elaboration_program_path = '~/ORAC/trunk/tools/primadorac/primadorac.bash',
                                                Protein_model = 0,
                                                Protein_chain = 'A',
                                                ph = 7.0,
                                                repairing_method = 'pdbfixer',
                                                MD_program = 'orac',
                                                MD_program_path = '~/ORAC/trunk/src/GNU-FFTW-OMP/orac',
                                                protein_tpg_file = '../amber99sb-ildn.tpg',
                                                protein_prm_file = '../amber99sb-ildn.prm',
                                                solvent_pdb = '../water.pdb',
                                                residue_substitution = 'standard',
                                                kind_of_processor = 'skylake',
                                                number_of_cores_per_node = 64,
                                                use_gpu = 'auto')

        if not os.path.exists(protein_id):
            os.mkdir(protein_id)
        os.chdir(protein_id)
        with open(f"{protein_id}_stdout.out", 'w') as stdout, contextlib.redirect_stdout(stdout):
            with open(f"{protein_id}_stderr.err", 'w') as stderr, contextlib.redirect_stderr(stderr):
                
                p.execute()

    except KeyboardInterrupt:
        break

    except Exception as err:
        print(err, err.args)

    else:
        print("Everything was ok")

    os.chdir('../.')

print("THIS...\nIS....\nTHE END")
