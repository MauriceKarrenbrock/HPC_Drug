import subprocess
import os
import HPC_Drug
import HPC_Drug.pipelines
import sys


protein_id = sys.argv[1]

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
                                        protein_tpg_file = 'amber99sb-ildn.tpg',
                                        protein_prm_file = 'amber99sb-ildn.prm',
                                        solvent_pdb = 'water.pdb')

# if not os.path.exists(protein_id):
#     os.mkdir(protein_id)
# os.chdir(protein_id)

p.execute()

#os.chdir('../.')
