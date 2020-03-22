"""
Copyright (c) 2020-2020 Maurice Karrenbrock

This software is open-source and is distributed under the
GNU Affero General Public License v3 (agpl v3) license

A copy of the license must be included with any copy of the program or part of it
"""

import subprocess

command_list = ['mol load pdb tests/2f3z.pdb\n',
                'package require autopsf\n',
                'autopsf << ENDMOL\n',
                'coordpdb /run/media/maurice/DATI/Documenti/maurice/università/Tesi_Magistrale/HPC_Drug/tests/2f3z.pdb\n',
                'patch\n'
                'guesscoord\n',
                'writepdb psf.pdb\n',
                'ENDMOL\n',
                'quit\n']

#pip = subprocess.Popen('vmd -dispdev text\n', shell = True)#, check = True)


# Split a file containing protein and water into separate segments.
# Creates files named myfile_water.pdb, myfile_frag0.pdb, myfile_frag1.pdb,...
# Requires VMD.
comandi = [
'mol load pdb myfile.pdb',
'set water [atomselect top water]',
'$water writepdb myfile_water.pdb',
'set protein [atomselect top protein]',
'set chains [lsort -unique [$protein get pfrag]]',
'foreach chain $chains {',
'set sel [atomselect top "pfrag $chain"]',
'$sel writepdb myfile_frag${chain}.pdb',
'}']

# substitutions = {
#     '2AS':'ASP', '3AH':'HIS', '5HP':'GLU', 'ACL':'ARG', 'AGM':'ARG', 'AIB':'ALA', 'ALM':'ALA', 'ALO':'THR', 'ALY':'LYS', 'ARM':'ARG',
#     'ASA':'ASP', 'ASB':'ASP', 'ASK':'ASP', 'ASL':'ASP', 'ASQ':'ASP', 'AYA':'ALA', 'BCS':'CYS', 'BHD':'ASP', 'BMT':'THR', 'BNN':'ALA',
#     'BUC':'CYS', 'BUG':'LEU', 'C5C':'CYS', 'C6C':'CYS', 'CAS':'CYS', 'CCS':'CYS', 'CEA':'CYS', 'CGU':'GLU', 'CHG':'ALA', 'CLE':'LEU', 'CME':'CYS',
#     'CSD':'ALA', 'CSO':'CYS', 'CSP':'CYS', 'CSS':'CYS', 'CSW':'CYS', 'CSX':'CYS', 'CXM':'MET', 'CY1':'CYS', 'CY3':'CYS', 'CYG':'CYS',
#     'CYM':'CYS', 'CYQ':'CYS', 'DAH':'PHE', 'DAL':'ALA', 'DAR':'ARG', 'DAS':'ASP', 'DCY':'CYS', 'DGL':'GLU', 'DGN':'GLN', 'DHA':'ALA',
#     'DHI':'HIS', 'DIL':'ILE', 'DIV':'VAL', 'DLE':'LEU', 'DLY':'LYS', 'DNP':'ALA', 'DPN':'PHE', 'DPR':'PRO', 'DSN':'SER', 'DSP':'ASP',
#     'DTH':'THR', 'DTR':'TRP', 'DTY':'TYR', 'DVA':'VAL', 'EFC':'CYS', 'FLA':'ALA', 'FME':'MET', 'GGL':'GLU', 'GL3':'GLY', 'GLZ':'GLY',
#     'GMA':'GLU', 'GSC':'GLY', 'HAC':'ALA', 'HAR':'ARG', 'HIC':'HIS', 'HIP':'HIS', 'HMR':'ARG', 'HPQ':'PHE', 'HTR':'TRP', 'HYP':'PRO',
#     'IAS':'ASP', 'IIL':'ILE', 'IYR':'TYR', 'KCX':'LYS', 'LLP':'LYS', 'LLY':'LYS', 'LTR':'TRP', 'LYM':'LYS', 'LYZ':'LYS', 'MAA':'ALA', 'MEN':'ASN',
#     'MHS':'HIS', 'MIS':'SER', 'MLE':'LEU', 'MPQ':'GLY', 'MSA':'GLY', 'MSE':'MET', 'MVA':'VAL', 'NEM':'HIS', 'NEP':'HIS', 'NLE':'LEU',
#     'NLN':'LEU', 'NLP':'LEU', 'NMC':'GLY', 'OAS':'SER', 'OCS':'CYS', 'OMT':'MET', 'PAQ':'TYR', 'PCA':'GLU', 'PEC':'CYS', 'PHI':'PHE',
#     'PHL':'PHE', 'PR3':'CYS', 'PRR':'ALA', 'PTR':'TYR', 'PYX':'CYS', 'SAC':'SER', 'SAR':'GLY', 'SCH':'CYS', 'SCS':'CYS', 'SCY':'CYS',
#     'SEL':'SER', 'SEP':'SER', 'SET':'SER', 'SHC':'CYS', 'SHR':'LYS', 'SMC':'CYS', 'SOC':'CYS', 'STY':'TYR', 'SVA':'SER', 'TIH':'ALA',
#     'TPL':'TRP', 'TPO':'THR', 'TPQ':'ALA', 'TRG':'LYS', 'TRO':'TRP', 'TYB':'TYR', 'TYI':'TYR', 'TYQ':'TYR', 'TYS':'TYR', 'TYY':'TYR'
# }
# proteinResidues = ['ALA', 'ASN', 'CYS', 'GLU', 'HIS', 'LEU', 'MET', 'PRO', 'THR', 'TYR', 'ARG', 'ASP', 'GLN', 'GLY', 'ILE', 'LYS', 'PHE', 'SER', 'TRP', 'VAL']
#rnaResidues = ['A', 'G', 'C', 'U', 'I']
#dnaResidues = ['DA', 'DG', 'DC', 'DT', 'DI']
protein_filename = 'tests/2f3z.pdb'
protein_id = '2f3z'

k = [
'# Split a file containing protein and water into separate segments.',
'# Creates files named myfile_water.pdb, myfile_frag0.pdb, myfile_frag1.pdb,...',
'# Requires VMD.',
f'mol load pdb {protein_filename}',
'set water [atomselect top water]',
f'$water writepdb {protein_id}_water.pdb',
'set protein [atomselect top protein]',
'set chains [lsort -unique [$protein get pfrag]]',
'foreach chain $chains {',
'       set sel [atomselect top "pfrag $chain"]',
'       $sel writepdb ' + protein_id + '_frag${chain}.pdb',
'}',
'',
'package require psfgen',
'foreach chain $chains {',
'       psfgen << ENDMOL',
'       topology toph19.inp',
'',
'       segment PROT_${chain} {',
'               pdb PROT${chain}_protein.pdb',
'       }',
'',
'       coordpdb PROT${chain}_protein.pdb PROT_${chain}',
'',
'       segment SOLV {',
'               auto none ',
'               pdb PROT${chain}_water.pdb',
'       }',
'',
'       coordpdb output/6PTI_water.pdb SOLV',
'',
'       guesscoord',
'',
'       writepsf PROT${chain}.psf',
'       writepdb PROT${chain}.pdb',
'',
'       ENDMOL',
'}'
]

# with open('psf_stript_py.tcl', 'w') as f:
#     for line in k:
#         f.write(line)
# import os
# pip = subprocess.run(f'vmd -dispdev text -e < {os.getcwd()}/psf_stript_py.tcl', shell = True, check = True)

# import pdbfixer
# import simtk.openmm.app
# fixer = pdbfixer.PDBFixer(filename='pdb2f3z.ent')
# fixer.findMissingResidues()
# fixer.findNonstandardResidues()
# fixer.replaceNonstandardResidues()
# fixer.removeHeterogens(False)
# fixer.findMissingAtoms()
# fixer.addMissingAtoms()
# fixer.addMissingHydrogens(7.0)
# #fixer.addSolvent(fixer.topology.getUnitCellDimensions())
# simtk.openmm.app.PDBFile.writeFile(fixer.topology, fixer.positions, open('Cazzo.pdb', 'w'))

# import Bio.PDB
# import Bio.PDB.MMCIF2Dict

# p = Bio.PDB.MMCIFParser()
# struct = p.get_structure('2rfh', '2rfh.cif')

# cif_dict = Bio.PDB.MMCIF2Dict.MMCIF2Dict('2rfh.cif')

# print(cif_dict['_entity_poly_seq.mon_id'])
# s = Bio.PDB.PDBIO()
# s.set_structure(struct)
# s.save('2rfh_bio.pdb')

import pdbfixer
import simtk.openmm.app

input_filename = '1df8.cif'
output_filename = 'Cazzo.pdb'

with open(input_filename, 'r') as f:
    fixer = pdbfixer.PDBFixer(pdbxfile = f)
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
#fixer.removeHeterogens(False)
fixer.findMissingAtoms()
fixer.addMissingAtoms()
#fixer.addMissingHydrogens(8.0)
#fixer.addSolvent(fixer.topology.getUnitCellDimensions())

with open(output_filename, 'w') as w:
    simtk.openmm.app.PDBFile.writeFile(fixer.topology, fixer.positions, open(output_filename, 'w'), keepIds= True)

# import importlib_resources
# import lib
# with importlib_resources.path('lib', 'amber99sb-ildn.tpg') as path:
#     print(path)
# print(path)


# def preprocess_source(self, in_file, additional_args=[]):
#         import subprocess

#         self._args.extend(self._build_compiler_flags())
#         self._args.extend(additional_args)

#         result = subprocess.run(self._args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

#         if result.returncode == 0:
#             return result.stdout
#         else:
#             if result.stderr:
#                 Style.error('Preprocess failed: ')
#                 print(result.stderr)

#             return '' 

#Seqres da header : _pdbx_poly_seq_scheme.mon_id