"""
This file contains the classes and methods to find structural informations making calculations
on the structure
like finding cys-cys bounds or the metal binding residues when there are no available
or not up to date
this methods are much slower than parsing a mmcif header
"""
from HPC_Drug import file_manipulation
from HPC_Drug import structures

import Bio.PDB


class StructureInformations(file_manipulation.FileCruncer):

    def __init__(self, Protein = None):
        
        self.Protein = Protein
        if not isinstance(self.Protein, structures.Protein):
            raise TypeError(f"Need a Protein instance not a {type(self.Protein)}")
    
    def get_structure(self, filename = None, file_type = "pdb", protein_id = "dummy"):
        """Parses a pdb or a mmcif with Biopython, returns a biopython structure"""

        if filename == None:
            filename = self.Protein.filename
            file_type = self.Protein.file_type
            protein_id = self.Protein.protein_id


        if file_type == "pdb":
            p = Bio.PDB.PDBParser()
        
        elif file_type == "cif":
            p = Bio.PDB.MMCIFParser()

        else:
            raise TypeError(f"file_type can only be 'pdb' or 'cif' types, not {file_type}")

        structure = p.get_structure(protein_id, filename)

        return structure
