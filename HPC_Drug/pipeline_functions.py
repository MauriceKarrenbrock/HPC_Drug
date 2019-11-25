from HPC_Drug import file_manipulation
from HPC_Drug import structures
from HPC_Drug import important_lists
import Bio.PDB
import Bio.PDB.MMCIF2Dict
import collections

def parse(Protein):
    """Parses the PDBx/mmCIF file
    returns a protein containing informations about
    the residues binding a metal, the organic ligand resname (if present)
    the sulfidic bonds and the seqres

    takes and resturns a structures.Protein instance and a list of
    ligand resnames and resnumbers [[resname, resnumber], [..., ...], ...]

    parses the seqres from the mmcif header if possible
    """

    #Parsing substitutions, sulf bonds and the resname of the organic ligand
    subst_parser = file_manipulation.SubstitutionParser()

    Protein.substitutions_dict, Protein.sulf_bonds, ligand_resnames =\
        subst_parser.parse_substitutions_PDB(file_name = Protein.filename,
                                            protein_chain = Protein.chain)
    
    ligand_resnames = subst_parser.get_ligand_resnum(Protein = Protein,
                                                    ligand_resnames = ligand_resnames,
                                                    chain_model_selection = True)

    try:
        Protein.seqres = get_seqres_mmcif_header(filename = Protein.filename)
    
    except:
        Protein.seqres = None
    
    return Protein, ligand_resnames


def custom_orac_seqres_from_PDB(Protein):
    """Gets the seqres from a pdb file
    and stores it in a protein instance

    param:: a Protein instance"""

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

def get_seqres_mmcif_header(Protein = None, filename = None):
    """Gets the seqres written in the mmCIF header
    can take a Protein instance or a filename

    returns a protein instance if givent a protein
    and a list if given a filename"""

    if Protein != None:
        cif_dict = Bio.PDB.MMCIF2Dict.MMCIF2Dict(Protein.filename)
        Protein.seqres = cif_dict['_pdbx_poly_seq_scheme.mon_id']
        return Protein

    elif filename != None:
        cif_dict = Bio.PDB.MMCIF2Dict.MMCIF2Dict(filename)
        seqres = cif_dict['_pdbx_poly_seq_scheme.mon_id']
        return seqres
    
    else:
        raise TypeError("Need a valid Protein or a valid filename")

def update_sulf_bonds(Protein = None):
    """Updates Protein.sulf_bonds using Protein.cys_dict
    takes a Protein instance and returns one"""

    #tmp_sulf_bonds contains the cysteine_numbers of the sulf_bonds
    tmp_sulf_bonds = []
    for bond in Protein.sulf_bonds:
        
        dummy_var = (Protein.cys_dict[bond[0]], Protein.cys_dict[bond[1]])
        tmp_sulf_bonds.append(dummy_var)

    if Protein.file_type == 'pdb':
        p = Bio.PDB.PDBParser()
    
    elif Protein.file_type == 'cid':
        p = Bio.PDB.MMCIFParser()

    #parsing the new pdb to check the new cys resnums
    struct = p.get_structure(Protein.protein_id, Protein.filename)

    residues = struct.get_residues()

    i = 0
    tmp_cys_dict = {}
    for residue in residues:
        if residue.resname in important_lists.cyst_resnames:
            
            #i = cysteine_number ; residue._id[1] = resnum
            i = i + 1
            tmp_cys_dict[str(i)] = str(residue._id[1])

    #Creating the new sulf_bonds
    new_sulf_bonds = []

    for bond in tmp_sulf_bonds:
        
        dummy_var = (tmp_cys_dict[bond[0]], tmp_cys_dict[bond[1]])
        new_sulf_bonds.append(dummy_var)

    Protein.sulf_bonds = new_sulf_bonds

    return Protein


def get_iterable(x):
    """Returns an iterable, even if given a single value"""
    if isinstance(x, collections.Iterable):
        return x
    else:
        return (x,)