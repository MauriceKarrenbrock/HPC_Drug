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
import important_lists
import pipeline_functions
import file_manipulation
import structures

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
        given the matrix containing the new basisi vectors as columns
        The default is to use the 
        base of the autovectors of the moment of inertia matrix of the structure
        
        new_vec = dot(M ** -1, old_vec)"""

        if structure == None:
            Pp = Bio.PDB.PDBParser()
            structure = Pp.get_structure(self.Protein.protein_id, self.Protein.filename)
        
        if rot_matrix == None:
            tmp, tmp_1, rot_matrix = self.calculate_moment_of_intertia_tensor()
        
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
        """Calculates the disatance between the center ofmass of two given structure
        returns the two centers of mass and their distance"""

        if structure_1 == None or structure_2 == None:
            raise TypeError('Need a Biopython strucure, not a None type')

        #Get the center of mass of both as COM objects
        COM_1 = self.center_of_mass(structure_1)
        COM_2 = self.center_of_mass(structure_2)

        distance = (COM_1.coord[0] - COM_2.coord[0])**2. + (COM_1.coord[1] - COM_2.coord[1])**2.
        distance = distance ** (0.5)

        #Return the centers of mass coordinates and the distance between them
        return COM_1.coord, COM_2.coord, distance
    
    def separate_protein_ligand(self, Protein = None, Ligand = None):
        """Takes a Protein instance whose filename is the pdbfile with both
        the protein and the ligands together and returns the separated biopython structures"""

        if Protein == None:
            Protein = self.Protein

        if Ligand == None:
            Ligand = self.Ligand

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
            ligand.structure = c.get_ligand(Protein.protein_id, Protein.filename, ligand.ligand_resname)
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
                if "ATOM" in line or "HETATM" in line:

                    line = line.split()

                    atom_number.append(line[1].strip())
                    resname.append(line[3].strip().upper())
        
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