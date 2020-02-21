# Contains the classes and functions to orient
#the given strucutre along the inertia's tensor
#Create a box around it and do many other calculations needed to create
#the input for the orac optimization with solvent

from Bio.PDB import Entity
import Bio.PDB
import Bio.PDB.vectors
from math import sqrt
from numpy import *
#import sys
from HPC_Drug import important_lists
from HPC_Drug import pipeline_functions
from HPC_Drug import file_manipulation
from HPC_Drug import structures

#deactivating all
#BiopythonWarning
import warnings
import Bio
warnings.simplefilter('ignore', Bio.BiopythonWarning)


class Orient(object):
    """This class contains the methods needed to 
    put the given stucture in a box with solvent
    and do many other calculations needed to create
    the input for the orac optimization with solvent"""

    def __init__(self, Protein = None, Ligand = None):
        self.Protein = Protein
        self.Ligand = Ligand
        self.atom_weights = important_lists.atom_weights


    def center_of_mass(self, entity = None, geometric=False):
        """
        Returns gravitic [default] or geometric center of mass of an Entity.
        Geometric assumes all masses are equal (geometric=True)
        """

        # Copyright (C) 2010, Joao Rodrigues (anaryin@gmail.com)
        # https://github.com/JoaoRodrigues/biopython/blob/GSOC2010/Bio/Struct/Geometry.py
        
        if entity == None:
            Pp = Bio.PDB.PDBParser()
            entity = Pp.get_structure(self.Protein.protein_id, self.Protein.filename)
        
        # Structure, Model, Chain, Residue
        if isinstance(entity, Entity.Entity):
            atom_list = entity.get_atoms()
        # List of Atoms
        elif hasattr(entity, '__iter__') and [x for x in entity if x.level == 'A']:
            atom_list = entity
        else: # Some other weirdo object
            raise ValueError("Center of Mass can only be calculated from the following objects:\n"
                                "Structure, Model, Chain, Residue, list of Atoms.")
        class COM:
            def __init__(self,coord):
                self.coord=coord

        positions = [ [], [], [] ] # [ [X1, X2, ..] , [Y1, Y2, ...] , [Z1, Z2, ...] ]
        masses = []

        for atom in atom_list:
            try:
                atom.mass = self.self.atom_weights[atom.name[0].capitalize()]
            except:
                atom.mass = 1.0

            masses.append(atom.mass)

            for i, coord in enumerate(atom.coord.tolist()):
                positions[i].append(coord)

        # If there is a single atom with undefined mass complain loudly.
        if 'ukn' in set(masses) and not geometric:
            raise ValueError("Some Atoms don't have an element assigned.\n"
                            "Try adding them manually or calculate the geometrical center of mass instead.")

        if geometric:
            com = COM([sum(coord_list)/len(masses) for coord_list in positions])
            return com
        else:
            w_pos = [ [], [], [] ]
            for atom_index, atom_mass in enumerate(masses):
                w_pos[0].append(positions[0][atom_index]*atom_mass)
                w_pos[1].append(positions[1][atom_index]*atom_mass)
                w_pos[2].append(positions[2][atom_index]*atom_mass)
            com = COM([sum(coord_list)/sum(masses) for coord_list in w_pos])
            return com

    def calculate_moment_of_intertia_tensor(self, structure = None):
        """
        Calculates the moment of inertia tensor from the molecule.
        Returns a numpy matrix.
        """

        #This is part of CCP4
        #http://www.ccp4.ac.uk/dist/checkout/arcimboldo/src/geometry.py

        if structure == None:
            Pp = Bio.PDB.PDBParser()
            structure = Pp.get_structure(self.Protein.protein_id, self.Protein.filename)

        if isinstance(structure, Entity.Entity):
            atom_list = structure.get_atoms()
            # List of Atoms
        elif hasattr(structure, '__iter__') and [x for x in structure if x.level == 'A']:
            atom_list = structure
        else: # Some other weirdo object
            raise ValueError("Center of Mass can only be calculated from the following objects:\n"
                                "Structure, Model, Chain, Residue, list of Atoms.")
        s_mass = 0.0
        for atom in atom_list:

            try:
                atom.mass = self.atom_weights[atom.name[0].capitalize()]
            except:
                atom.mass = self.atom_weights[atom.element.capitalize()]
            s_mass += atom.mass

        com = self.center_of_mass(structure, False)
        cx, cy, cz = com.coord

        n_atoms = 0
        tensor_xx, tensor_xy, tensor_xz = 0, 0, 0
        tensor_yx, tensor_yy, tensor_yz = 0, 0, 0
        tensor_zx, tensor_zy, tensor_zz = 0, 0, 0
        #s_mass = sum([a.mass for a in atom_list])

        if isinstance(structure, Entity.Entity):
            atom_list = structure.get_atoms()
        elif hasattr(structure, '__iter__') and [x for x in structure if x.level == 'A']:
            atom_list = structure

        for atom in atom_list:
            ax, ay, az = atom.coord
            try:
                tensor_xx += ((ay-cy)**2 + (az-cz)**2)*self.atom_weights[atom.name[0].capitalize()]
                tensor_xy += -1*(ax-cx)*(ay-cy)*self.atom_weights[atom.name[0].capitalize()]
                tensor_xz += -1*(ax-cx)*(az-cz)*self.atom_weights[atom.name[0].capitalize()]
                tensor_yy += ((ax-cx)**2 + (az-cz)**2)*self.atom_weights[atom.name[0].capitalize()]
                tensor_yz += -1*(ay-cy)*(az-cz)*self.atom_weights[atom.name[0].capitalize()]
                tensor_zz += ((ax-cx)**2 + (ay-cy)**2)*self.atom_weights[atom.name[0].capitalize()]
            except:
                tensor_xx += ((ay-cy)**2 + (az-cz)**2)*self.atom_weights[atom.element.capitalize()]
                tensor_xy += -1*(ax-cx)*(ay-cy)*self.atom_weights[atom.element.capitalize()]
                tensor_xz += -1*(ax-cx)*(az-cz)*self.atom_weights[atom.element.capitalize()]
                tensor_yy += ((ax-cx)**2 + (az-cz)**2)*self.atom_weights[atom.element.capitalize()]
                tensor_yz += -1*(ay-cy)*(az-cz)*self.atom_weights[atom.element.capitalize()]
                tensor_zz += ((ax-cx)**2 + (ay-cy)**2)*self.atom_weights[atom.element.capitalize()]


        in_tensor =  array([[tensor_xx, tensor_xy, tensor_xz], [tensor_xy, tensor_yy, tensor_yz], [tensor_xz,
            tensor_yz, tensor_zz]])
        D,V = linalg.eig(in_tensor)

        a = sqrt((5/(2*s_mass)) * (D[0] - D[1] + D[2]))
        b = sqrt((5/(2*s_mass)) * (D[2] - D[0] + D[1]))
        c = sqrt((5/(2*s_mass)) * (D[1] - D[2] + D[0]))
        return sorted([a, b, c]),D,V

    def base_change_structure(self, structure = None, rot_matrix = None):
        """Changes the base of the coordinates of a given Biopython strucure
        given the matrix containing the new basis vectors as columns
        The default is to use the 
        base of the autovectors of the moment of inertia matrix of the structure
        
        new_vec = dot(M ** -1, old_vec)"""

        if structure == None:
            Pp = Bio.PDB.PDBParser()
            structure = Pp.get_structure(self.Protein.protein_id, self.Protein.filename)
        
        if rot_matrix == None:
            tmp, tmp_1, rot_matrix = self.calculate_moment_of_intertia_tensor()

        #I want to be sure that the matrix doesn't contain a reflection (det < 0)
        #because it would messup the chirality
        #In case there is one i reflect the matrix around the z axis
        if linalg.det(rot_matrix) < 0:
            rot_matrix[:,:3] = - rot_matrix[:,:3]
        
        #make the inverse
        rot_matrix = linalg.inv(rot_matrix)
        
        atoms = structure.get_atoms()
    
        for atom in atoms:
            
            #Rotating the vector in the new basis
            atom.coord = rot_matrix @ atom.coord


        return structure

    def create_box(self, rotated_structure = None, border = 16.0):
        """Creates a box that contains the given structure
        The structure shall be rotated along the moment of inertia matrix autovectors"""

        Max_x, min_x = None, None
        Max_y, min_y = None, None
        Max_z, min_z = None, None

        atoms = rotated_structure.get_atoms()

        #Giving a starting value t max and min
        Max_x, Max_y, Max_z = -1.E-10, -1.E-10, -1.E-10
        min_x, min_y, min_z = 1.E10, 1.E10, 1.E10

        for atom in atoms:
            
            x, y, z = atom.coord
            
            Max_x = max(Max_x, x)
            Max_y = max(Max_y, y)
            Max_z = max(Max_z, z)

            min_x = min(min_x, x)
            min_y = min(min_y, y)
            min_z = min(min_z, z)

        #The box will be border/2 larger then the protein along any side
        lx = Max_x - min_x + border
        ly = Max_y - min_y + border
        lz = Max_z - min_z + border

        return lx, ly, lz

    def create_solvent_grid(self, lx = None, ly = None, lz = None, solv_dist = 3.6):
        """Generates the number of waters to put in the box
        of dimentions lx, ly, lz
        with solvent distance = solv dist (default = 3.6 A)
        
        returns the number of solvent molucules along each axis nx, ny, nz
        (int)"""

        nx = int(lx / solv_dist)

        ny = int(ly / solv_dist)

        nz = int(lz / solv_dist)

        return nx, ny, nz
    
    def create_recipr_solvent_grid(self, lx = None, ly = None, lz = None, solv_dist = 3.6):
        """Cretes the informations about the grid on the reciprocal
        reticle, needed for the FFT on atomic charges
        
        EWALD PME"""

        pme_x, pme_y, pme_z = self.create_solvent_grid(lx, ly, lz, solv_dist)


        pme_x *= 3
        pme_y *= 3
        pme_z *= 3

        pme = [pme_x, pme_y, pme_z]

        #I find the nearest multiple of 2, 5, 7
        for i in range(0, len(pme)):
                while True:
                    if pme[i] % 2 == 0:
                        break
                    elif pme[i] % 5 == 0:
                        break
                    elif pme[i] % 7 == 0:
                        break
                    else:
                        pme[i] = pme[i]-1
        return pme

    def center_mass_distance(self, structure_1 = None, structure_2 = None):
        """Calculates the disatance between the center of mass of two given structures
        returns the two centers of mass and their distance"""

        if structure_1 == None or structure_2 == None:
            raise TypeError('Need a Biopython strucure, not a None type')

        #Get the center of mass of both as COM objects
        COM_1 = self.center_of_mass(structure_1)
        COM_2 = self.center_of_mass(structure_2)

        distance = (COM_1.coord[0] - COM_2.coord[0])**2. + (COM_1.coord[1] - COM_2.coord[1])**2. + (COM_1.coord[2] - COM_2.coord[2])**2.
        distance = distance ** (0.5)

        #Return the centers of mass coordinates and the distance between them
        return COM_1.coord, COM_2.coord, distance
    
    def separate_protein_ligand(self, Protein = None, Ligand = None):
        """Takes a Protein instance whose filename is the pdbfile with both
        the protein and the ligands together and returns the separated biopython structures"""
        
        def update_ligand_resnumbers(Protein = self.Protein, Ligand = self.Ligand):
            """Updates the ligand resnumbers that may be changed
            because of previous calculations"""

            if Ligand == None or len(Ligand) == 0:
                print("I found no ligand, returning None")
                return None

            ligand_resnames = []
            for ligand in pipeline_functions.get_iterable(Ligand):
                if ligand.ligand_resname not in ligand_resnames:
                    ligand_resnames.append(ligand.ligand_resname)

            c = file_manipulation.SubstitutionParser()
            ligand_resnames = c.get_ligand_resnum(Protein = Protein,
                                                ligand_resnames = ligand_resnames,
                                                chain_model_selection = False)
            #now ligand resnames is in the form [[ligand_resname, ligand_resnumber], [..., ...], ...]

            #Creating a completely new Ligand list
            Ligand = []

            for ligand in ligand_resnames:
                Ligand.append(structures.Ligand(ligand_resname = ligand[0],
                                                filename = None,
                                                structure = None,
                                                file_type = 'pdb',
                                                topology_file = None,
                                                param_file = None,
                                                res_number = ligand[1]))

            return Ligand


        if Protein == None:
            Protein = self.Protein

        if Ligand == None:
            Ligand = self.Ligand

        Ligand = update_ligand_resnumbers(Protein = Protein, Ligand = Ligand)

        #Extract the protein with ProDy, write a new pdb
        #and then get the Biopython protein structure
        c = file_manipulation.PDBCruncer()
        tmp_protein = structures.Protein(protein_id= Protein.protein_id)
        tmp_protein.structure = c.get_protein(Protein.protein_id, Protein.filename)
        tmp_protein.filename = tmp_protein.write_PDB(f"{Protein.protein_id}_protein.pdb")

        p = Bio.PDB.PDBParser()
        
        Protein.structure = p.get_structure(tmp_protein.protein_id, tmp_protein.filename)

        #doing the same thing for every ligand
        ligand_structures = []
        for i, ligand in enumerate(pipeline_functions.get_iterable(Ligand)):
            ligand.structure = c.get_ligand(protein_id = Protein.protein_id,
                                            filename = Protein.filename,
                                            ligand_name =ligand.ligand_resname,
                                            ligand_resnumber = ligand.res_number)

            ligand.filename = ligand.write_PDB(f"{Protein.protein_id}_lgand{i}.pdb")
            ligand.structure = p.get_structure(Protein.protein_id, ligand.filename)
            ligand_structures.append(ligand.structure)

        return Protein.structure, ligand_structures
    
    def atom_numbers(self, Protein = None, Ligand = None):
        """Given the instance of a ligand and of a protein containing
        the filename of a pdb containing protein + ligands
        will give the beginning and end atom number of the protein
        and of any organic ligand
        
        return protein_atom, ligand atom
        protein_atom = [MAX, min]
        ligand_atom = [[MAX, min], [MAX, min], ...]"""

        if Protein == None:
            Protein = self.Protein
        
        if Ligand == None:
            Ligand = self.Ligand
        
        #Reading the PDB and getting the atom numbers and resnames of any atom
        atom_number = []
        resname = []

        with open(Protein.filename, 'r') as f:
            for line in f:
                if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':

                    atom_number.append(line[4:11].strip())
                    resname.append(line[17:20].strip().upper())
        
        #getting a list of the resnames of the ligands
        ligand_resnames = []

        for ligand in pipeline_functions.get_iterable(Ligand):
            ligand_resnames.append(ligand.ligand_resname)
        
        ligand_atoms = []

        protein_atoms = []
        
        #Maximum and minimum atom number for the protein
        M_p = int(-1E-10)
        m_p = int(1E10)

        for n, i in enumerate(atom_number):

            if resname[n] not in ligand_resnames:
                M_p = max(int(i), M_p)
                m_p = min(int(i), m_p)
        
        protein_atoms = [M_p, m_p]


        #Maximum and minimum atom number for any given ligand
        M_l = int(-1E-10)
        m_l = int(1E10)


        for ligand in pipeline_functions.get_iterable(Ligand):

            for n, i in enumerate(atom_number):

                if resname[n] == ligand.ligand_resname:
                    M_l = max(int(i), M_l)
                    m_l = min(int(i), m_l)
            
            ligand_atoms.append([M_l, m_l])

        return protein_atoms, ligand_atoms

    def get_hot_residues_for_rem(self, Protein = None, Ligand = None, cutoff = 3.0, residue_dist = 8.0):

        """ This functions will give the residue_id of any "hot residue" needed for the REM
        symulations, the hot residues are the one in contact with the given organic ligand

        It takes a Protein instance and a Ligand instance (or a list of them) but if the list
        is longer than one it will raise an error

        The protein instance must contain a file with both the protein and the ligand

        It returns a list of the redisue_id of the right residues

        cutoff is a float that tells the cutoff distance between the ligand and the nearest residue atom
        to be considered in contact default = 3.0 angstrom

        residue dist tells the maximum distance that a ligand atom and a residue center of mass
        must have to be scanned for the cutoff default 8.0 angstrom

        returns a list of listst: [[residue_resname, residue_id], ...]
        """

        if Protein == None:
            Protein = self.Protein
        
        if Ligand == None:
            Ligand = self.Ligand

        #as Ligand is often a list I check that there is only one
        Ligand = pipeline_functions.get_iterable(Ligand)
        if len(Ligand) != 1:
            raise ValueError(f"Can process only a Ligand at the time not {len(Ligand)}")
        Ligand = Ligand[0]

        if Protein.file_type == 'cif': 
            p = Bio.PDB.MMCIFParser()

        elif Protein.file_type == 'pdb':
            p = Bio.PDB.PDBParser()

        else:
            raise TypeError(f'"{Protein.file_type}" is not a valid type\n only "cif" and "pdb"')

        structure = p.get_structure(Protein.protein_id, Protein.filename)
        structure = structure[Protein.model][Protein.chain]

        #I refresh the ligand's residue id
        Sub_Parser = file_manipulation.SubstitutionParser()
        TMP_list = Sub_Parser.get_ligand_resnum(Protein = Protein,
                                                    ligand_resnames = Ligand.ligand_resname,
                                                    chain_model_selection = True)
        Ligand.res_number = TMP_list[0][1]

        output_list = []

        #for any atom of the ligand I iterate through the whole protein and first check for the nearest
        #residues then for the nearest atom of the residue, if it's distance is <= cutoff
        #it is a hot residue
        for atom in structure[Ligand.res_number]:
            
            for residue in structure:
                
                #checking if it's a valid residue
                is_valid_residue = residue.id[1] != Ligand.res_number and residue.resname not in important_lists.metals and residue.resname not in important_lists.trash and residue.id[0].strip() == ''
                if is_valid_residue:

                    COM_ligand_atom, COM_residue, distance = self.center_mass_distance(atom, residue)

                    if distance <= residue_dist:

                        TMP_atom_dist = 1.E+20
                        #check for the nearest atom of the residue
                        for other_atom in residue:

                            d = (atom.coord[0] - other_atom.coord[0])**2. + (atom.coord[1] - other_atom.coord[1])**2. + (atom.coord[2] - other_atom.coord[2])**2.
                            d = d ** (0.5)
                            
                            if d < TMP_atom_dist:
                                
                                TMP_atom_dist = d
                        
                        #checking if the nearest atom is near enough to be part of a hot residue
                        if TMP_atom_dist <= cutoff:
                            #I append the residue resname, the nearest atom name, and the residue id to the output list
                            output_list.append([residue.resname.strip().upper(), residue.id[1]])

        #useless variables
        COM_ligand_atom = None
        COM_residue = None


        return output_list