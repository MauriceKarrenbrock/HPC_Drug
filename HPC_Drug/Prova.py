import structures
import file_manipulation as fm
import os

filename = fm.download_protein_structure('2f3z', 'pdb', os.getcwd())

cruncer = fm.ProteinCruncer('pdb')

Protein = cruncer.get_protein('2f3z', filename, None)

Ligand = cruncer.get_ligand('2f3z', filename, '2f3z_ligand', None)

Protein.write_PDB(filename + '.prova')
Ligand.write_PDB('2f3z_ligand.pdb')

