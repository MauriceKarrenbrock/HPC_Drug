from HPC_Drug import file_manipulation
from HPC_Drug import structures
import collections

def parse(Protein):
    """Parses the PDBx/mmCIF file
    returns a protein containing informations about
    the residues binding a metal, the organic ligand resname (if present)
    the sulfidic bonds and the seqres

    takes and resturns a structures.Protein instance and a Ligand one
    
    the protein mmcif is converted in a pdb, and Protein.filename is updated"""

    #Parsing substitutions, sulf bonds and the resname of the organic ligand
    subst_parser = file_manipulation.SubstitutionParser()

    Protein.substitutions_dict, Protein.sulf_bonds, ligand_resnames =\
        subst_parser.parse_substitutions_PDB(file_name = Protein.filename,
                                            protein_chain = Protein.chain)
    
    return Protein, ligand_resnames


def get_seqres_PDB(Protein):
    """Gets the seqres from a pdb file
    and stores it in a protein instance

    param:: a Protein instance"""

    import Bio.PDB

    p = Bio.PDB.PDBParser()
    struct = p.get_structure(Protein.protein_id, Protein.filename)

    Protein.seqres = []
    tmp_metals = []

    for model in struct:
        for chain in model:

            #get residue list of every chain
            tmp_seqres = []
            
            res_list = Bio.PDB.Selection.unfold_entities(chain, 'R')

            for residue in res_list:
                
                #skipping HETATMs
                if residue._id[0].strip() == '':
                    #In the pdb this residues are called CYM because Biopython structures can
                    #only have 3 letters resnames, but Orac's tpg file calls it CYSM
                    #So I am modifying it in the seqres list
                    #It's probably not the best way to do it
                    if residue.resname.strip() == 'CYM':
                        tmp_seqres.append('CYSM')
                    else:
                        tmp_seqres.append(residue.resname.strip())
                
                #putting metal ions and hetatms in a separated list
                else:
                    tmp_metals.append(residue.resname.strip())

            
            #rename first and last residue of every chain
            tmp_seqres[0] = tmp_seqres[0] + '-H'
            tmp_seqres[-1] = tmp_seqres[-1] + '-O'

            for tmp in tmp_seqres:
                Protein.seqres.append(tmp)

        for i, resname in enumerate(Protein.seqres):

            if resname == 'CYS':
                if not (str(i+1) in Protein.substitutions_dict.keys()):

                    Protein.seqres[i] = 'CYSH'
    
    for tmp in tmp_metals:
        Protein.seqres.append(tmp)
    
    return Protein

def get_iterable(x):
    """Returns an iterable, even if given a single value"""
    if isinstance(x, collections.Iterable):
        return x
    else:
        return (x,)