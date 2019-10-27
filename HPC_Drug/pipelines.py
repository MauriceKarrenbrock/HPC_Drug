# This file contains the possible pipelines depending from input
class Pipeline(object):
    # Gives Not Implemented error
    def __init__(self):
        raise Exception(NotImplementedError)

class ChoosePipeline(Pipeline):
    """Takes the pipeline input
    and selects the right pipeline to follow"""
    
    def __init__(self):
        pass

    def choose_pipeline(self, protein, protein_filetype = 'cif', ligand = None):
        pass

class CifNoLigand_Pipeline(Pipeline):
    """Protein is given as a PDBx/mmCIF
    The ligand is not given as input
    Will be searched in the protein file"""
    pass

class CifLigand_Pipeline(Pipeline):
    """Protein is given as a PDBx/mmCIF
    The ligand is given as input"""
    pass


class PdbNoLigand_Pipeline(Pipeline):
    """Protein is given as a PDB
    The ligand is not given as input
    Will be searched in the protein file"""
    pass

class PdbLigand_Pipeline(Pipeline):
    """Protein is given as a PDB
    The ligand is given as input"""
    pass