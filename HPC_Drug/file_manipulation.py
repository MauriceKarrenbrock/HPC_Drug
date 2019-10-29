# Contains the classes for file manipulation like PDB & mmCIF

def download_protein_structure(protein_id, file_type = 'pdb', pdir = None):
    """The function downloads a PDB or a mmCIF from wwwPDB in a selected directory
    the default directory is the working directory
    it returns the filename (str)"""
    import Bio.PDB
    import os

    if file_type == 'cif':
        _file_type = 'mmCif'
    else:
        _file_type = file_type
    
    if pdir == None:
        pdir = os.getcwd()

    pdbl = Bio.PDB.PDBList()
    filename = pdbl.retrieve_pdb_file(protein_id, False, pdir, _file_type, False)
    
    if not os.path.exists(filename):
        raise FileNotFoundError(f'Was not able to download the protein or to find {filename}')
    
    else:
        return filename


class FileCruncer(object):
    """Generic file cruncer class"""
    def __init__(self):
        raise NotImplementedError('This kind file was not implemented yet')

def ProteinCruncer(file_type):
    """Returns the right file cruncer class
    for the given file_type"""

    if file_type == 'pdb':
        return PDBCruncer()
    elif file_type == 'cif':
        return MMCIFCruncer()
    else:
        return FileCruncer()


class PDBCruncer(FileCruncer):
    """Contains the metods for working on PDB
    with the ProDy module"""

    def __init__(self):
        pass

    def parse(self, protein_id = None, filename = None, PDB_model = None):
        """Parses a PDB with the ProDy parser"""
        import prody
        
        k = {'model' : PDB_model}
        if filename == None:
            parser = prody.parsePDB(protein_id, **k)
        else:
            
            parser = prody.parsePDB(filename, **k)

        return parser
        
    def get_protein(self, protein_id = None, filename = None, PDB_model = None, structure = None):
        import structures
        import prody

        if structure == None:
            parser = self.parse(protein_id, filename, PDB_model)
        else:
            parser = structure

        protein = structures.Protein(protein_id, filename, None, 'pdb', PDB_model)

        protein.structure = parser.select('protein')

        return protein

    def get_ligand(self, protein_id = None, filename = None, ligand_name = None, PDB_model = None, structure = None):
        """Creates a PDB with all the not water HETATM
        with ProDy"""
        import structures
        import prody
        if structure == None:
            parser = self.parse(protein_id, filename, PDB_model)
        else:
            parser = structure

        ligand = structures.Ligand(ligand_name, filename, None, 'pdb')

        ligand.structure = parser.select('not protein and not water')

        return ligand
    
    def write_PDB(self, structure, pdb_name):
        """
        Write a prody protein to a pdb file
        :param protein: protein object from prody
        :param pdb_name: base name for the pdb file
        :return: None
        """
        import prody
        
        prody.writePDB(f"{pdb_name}", structure)



class MMCIFCruncer(FileCruncer):
    """Contains the metods for working on mmCIF
    with th ProDy module"""
    import prody

    def __init__(self):
        pass

    def parse(self, protein_id = None, filename = None, PDB_model = None):
        """Parses a PDB with the ProDy parser"""
        import prody
        
        if filename == None:
            parser = prody.parseCIF(protein_id, model = PDB_model)
        else:
            parser = prody.parseCIF(filename, model = PDB_model)

        return parser
    
    def get_protein(self, protein_id = None, filename = None, PDB_model = None, structure = None):
        import structures
        import prody

        if structure == None:
            parser = self.parse(protein_id, filename, PDB_model)
        else:
            parser = structure

        p = PDBCruncer()

        protein = p.get_protein(protein_id, filename, PDB_model, parser)

        return protein

    def get_ligand(self, protein_id = None, filename = None, ligand_name = None, PDB_model = None, structure = None):
        """Creates a PDB with all the not water HETATM
        with ProDy"""
        import structures
        import prody
        if structure == None:
            parser = self.parse(protein_id, filename, PDB_model)
        else:
            parser = structure

        p = PDBCruncer()

        ligand = p.get_ligand(protein_id, filename, PDB_model, parser)

        return ligand
    
    def write_MMCIF(self, structure, pdb_name):
        """
        Not implementes, writes a PDB instead
        using PDBCruncer() class
        """
        print("Not implemented, printing PDB instead")

        p = PDBCruncer()
        p.write_PDB(f"{pdb_name}", structure)
