######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################


import HPC_Drug.file_manipulation
import HPC_Drug.important_lists
import sys
import os
import Bio.PDB

p = Bio.PDB.MMCIFParser()
with open("not_downloaded.txt", 'a') as err:
    err.write("Not downloaded proteins\n")

with open("possible_trash.txt", 'w') as w:
    for line in sys.stdin:
        protein_id = line.split()[0].strip()
        ligname = line.split()[0].replace('(', '').replace(')', '').strip()

        try:
        
            filename = HPC_Drug.file_manipulation.download_protein_structure(protein_id, 'cif')
        except KeyboardInterrupt:
            break
        except:
            with open("not_downloaded.txt", 'a') as err:
                err.write(f"{protein_id}\n")
            continue

        struct = p.get_structure(protein_id, filename)
        residues = Bio.PDB.Selection.unfold_entities(struct, 'R')
        already_found = []

        for residue in residues:
            if residue.id[0].strip() != '':
                if residue.resname.strip() != ligname:
                    if residue.resname.strip() not in HPC_Drug.important_lists.metals:
                        if residue.resname.strip() not in HPC_Drug.important_lists.trash:
                            if residue.resname.strip() != "HOH":
                                if residue.resname.strip() not in already_found:
                                    already_found.append(residue.resname.strip())
                                    w.write(f'"{residue.resname.strip()}",\n')

        os.system(f"rm {filename}")
    
print("END")





