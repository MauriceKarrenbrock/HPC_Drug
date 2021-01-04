######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains the super class of the Orac input templates
"""

from HPC_Drug.structures import ligand
from HPC_Drug.structures import protein
from HPC_Drug.structures import update_ligands
from HPC_Drug import orient
from HPC_Drug import important_lists
from HPC_Drug import lib
from HPC_Drug.auxiliary_functions import get_iterable
from HPC_Drug.auxiliary_functions import path
from HPC_Drug.files_IO import write_on_files
from HPC_Drug.files_IO import read_file
from HPC_Drug.PDB import add_chain_id
from HPC_Drug.PDB import prody
from HPC_Drug.PDB import biopython
from HPC_Drug.PDB.structural_information import mmcif_header
from HPC_Drug.MD.orac import get_resnum_cutoff


import os
import importlib_resources
import subprocess

class OracInput(object):
    """
    Super class of any Orac input template object
    
    The only public methos are the constructor and execute()
    execute returns a HPC_Drug.structires.protein.Protein instance 
    """

    def __init__(self,
                Protein,
                solvent_pdb = None,
                MD_program_path = 'orac'):
        
        """
        Protein :: HPC_Drug.structures.protein.Protein instance

        solvent_pdb :: string, it is the pdb file that contains the coordinates of a solvent molecule 
        it is needed if there has to be added a solvent box around the protein 
        default HPC_Drug.lib "water.pdb"

        MD_program_path :: string, the absolute path to the orac executable 
        dafault will look for an executable called orac in the PATH and the working directory (in this order)
        """
        
        self.Protein = Protein
        
        self.orac_in_file = os.getcwd() + f"/{self.Protein.protein_id}_orac.in"

        self.solvent_pdb = solvent_pdb
        #if no path is given searches the standard water.pdb inside lib module
        if self.solvent_pdb == None:
            with importlib_resources.path('HPC_Drug.lib', 'water.pdb') as _path:
                self.solvent_pdb = str(_path.resolve())
        
        self.MD_program_path = path.absolute_programpath(program = MD_program_path)

        self.output_pdb_file = os.getcwd() + f"/{self.Protein.protein_id}_orac.pdb"

        #an instance of orient.Orient class
        self.orient = orient.Orient(self.Protein, self.Protein.get_ligand_list())

        self.template = []

        #some values that are needed many times are calculated in the constructor

        #The box sizes (lx, ly, lz)
        #the structure is rotated in its tensor of inertia ref
        self.box = self._create_selfbox()

    def _create_selfbox(self):
        """
        private

        This method is called by the constructor
        """

        #The box sizes (lx, ly, lz)
        #the structure is rotated in its tensor of inertia ref
        self.Protein.update_structure("biopython")
        self.Protein.structure = self.orient.base_change_structure()
        self.Protein.write(file_name = self.Protein.pdb_file, struct_type = 'biopython')
        box = self.orient.create_box(self.Protein.structure)

        return box

    def _write_box(self):
        """
        private

        Writes the string about the box size
        and rotates the strucure in a smart way
        """
        
        lx, ly, lz = self.box

        string = f"   CRYSTAL    {lx:.2f}    {ly:.2f}    {lz:.2f}  !! simulation BOX"

        return string

    def _write_sulf_bond_string(self):
        """
        private
        
        Writes the part of the template inherent to the disulf bonds
        """

        Protein = self.Protein

        cutoff = self._get_protein_resnumber_cutoff()
        
        tmp_sulf = []
        
        for bond in Protein.sulf_bonds:

            #Correcting with the right cutoff
            bound_1 = bond[0] + cutoff
            bound_2 = bond[1] + cutoff

            string = f"   bond 1sg 2sg residue   {bound_1}  {bound_2}" + "\n"
            #needs redundancy
            string = string + f"   bond 1sg 2sg residue   {bound_2}  {bound_1}"
            tmp_sulf.append(string)
        
        string = '\n'.join(tmp_sulf)

        return string

    def _get_protein_resnumber_cutoff(self):
        """
        private
        
        Orac starts the residue count from one
        but pdb files often don't
        so I get the id of the resnumber of the first residue
        
        returns the cutoff to apply to start from one
        """

        self.Protein.update_structure(struct_type = "biopython")

        cutoff = get_resnum_cutoff.get_resnum_cutoff(structure = self.Protein.structure)

        return cutoff

    def _write_ligand_tpg_path(self):
        """
        private
        
        Writes the tpg file for any given ligand
        """

        Ligand = self.Protein.get_ligand_list()

        if Ligand == [] or Ligand == None:
            return ""

        tmp = []
        for ligand in Ligand:

            string = f"   READ_TPG_ASCII {ligand.tpg_file} !! ligand"

            tmp.append(string)

        string = '\n'.join(tmp)

        return string

    def _write_ligand_prm_path(self):
        """
        private
        
        Writes the prm file for any given ligand
        """

        Ligand = self.Protein.get_ligand_list()

        if Ligand == [] or Ligand == None:
            return ""

        tmp = []
        for ligand in Ligand:

            string = f"   READ_PRM_ASCII {ligand.prm_file}  !! ligand"

            tmp.append(string)

        string = '\n'.join(tmp)

        return string

    def _get_ligand_name_from_tpg(self):
        """
        private
        
        Gets the name of the ligand from the tpg file
        for any given ligand
        """

        Ligand = self.Protein.get_ligand_list()

        if Ligand == [] or Ligand == None:
            return ""


        residue_strings = []

        for ligand in Ligand:
            
            lines = read_file.read_file(file_name = ligand.tpg_file)
            
            for line in lines:

                if 'RESIDUE' in line:
                    
                    string = f"      {line.split()[1].strip()} !! ligand name in tpg file"

                    residue_strings.append(string)

                    break

            else:
                raise RuntimeError(f'Could not find the residue name in {ligand.tpg_file}') 

        return '\n'.join(residue_strings)

    def _write_solvent_grid(self):
        """
        private
        
        Writes the informations about the solvent grid
        """

        lx, ly, lz = self.box

        nx, ny, nz = self.orient.create_solvent_grid(lx, ly, lz)

        string = f"   GENERATE   {nx}  {ny}  {nz}  !! the grid depends on the BOX"

        return string

    def _write_EWALD_PME(self):
        """
        private
        
        Writes the informations about the grid
        in the reciprocal reticle
        """
        
        lx, ly, lz = self.box

        pme_x, pme_y, pme_z = self.orient.create_recipr_solvent_grid(lx, ly, lz)

        string = f"   EWALD PME 0.37   {pme_x}  {pme_y}  {pme_z}   4 !! grid on the reciprocal"

        return string

    def _write_ADD_STR_COM(self):
        """
        private
        
        Writes the informations about the bounding between the protein and the ligands
        """

        Ligand = self.Protein.get_ligand_list()

        if Ligand == [] or Ligand == None:
            return ""


        #updates the Protein._ligands
        self.Protein = update_ligands.update_ligands(Protein = self.Protein)

        #get a list of ligand structures
        ligand_structures = []

        Ligand = self.Protein.get_ligand_list()
        for lig in Ligand:

            lig.update_structure(struct_type = "biopython")

            ligand_structures.append(lig.structure)

        self.Protein.update_structure(struct_type = "prody")

        #get the only protein biopython structure
        prody_select_obj = prody.ProdySelect(structure = self.Protein.structure)
        protein_structure = prody_select_obj.protein_and_ions()
        prody.write_pdb(structure = protein_structure, file_name = f"{self.Protein.protein_id}_protein.pdb")
        protein_structure = biopython.parse_pdb(protein_id = self.Protein.protein_id, file_name = f"{self.Protein.protein_id}_protein.pdb")

        
        #Get the atom number range of each segment
        protein_atoms, ligand_atoms = self.orient.protein_ligand_atom_numbers()

        string = []

        #get the distance of the centers of mass of the protein and the ligands
        #and create the string
        for i, ligand in enumerate(get_iterable.get_iterable(ligand_structures)):

            COM_Protein, COM_ligand, distance = self.orient.center_mass_distance(protein_structure, ligand)
            
            #These are useless
            COM_Protein = COM_ligand = None

            tmp_string = ["   ADD_STR_COM   !! linker COM-COM legand protein",
                        f"       ligand     {ligand_atoms[i][1]}      {ligand_atoms[i][0]}",
                        f"       target    {protein_atoms[1]}      {protein_atoms[0]}",
                        f"       force   0.15     {distance}    {distance+0.0001}",
                        "   END"]
            
            tmp_string = '\n'.join(tmp_string)
            string.append(tmp_string)

        if len(string) > 1:
            string = '\n'.join(string)
        elif len(string) == 0:
            string = ''
        else:
            string = string[0]
        
        return string

    def _write_LINKED_CELL(self):
        """
        private

        writes the LINKED CELL string
        """

        lx, ly, lz = self.box

        string = f"   LINKED_CELL   {int(lx/3.0)}  {int(ly/3.0)}  {int(lz/3.0)}"

        return string

    def _write_template_on_file(self, template = None):
        """
        private

        Writes the objects template on {filename} file
        """

        if template is None:
            template = self.template

        lines = ["\n".join(template)]
        
        write_on_files.write_file(lines = lines, file_name = self.orac_in_file)


    def _write_chain_in_pdb(self):
        """
        private

        This is a patch because orac removes the chain id from
        pdb files and this confuses some pdb parsers
        """

        add_chain_id.add_chain_id(pdb_file = self.output_pdb_file, chain = self.Protein.chain)

    def _run(self):
        """
        private
        """

        print("Running Orac")

        #Orac reads the input from the stdin
        with open(self.orac_in_file, 'r') as input_file:
            r = subprocess.run([self.MD_program_path],
                            shell = False,
                            stdout = subprocess.PIPE,
                            stderr = subprocess.PIPE,
                            universal_newlines = True,
                            stdin = input_file)

        print(r.stdout)
        print(r.stderr)

        if r.returncode != 0:
            raise Exception(f"Orac failure\n{r.stdout}\n{r.stderr}")

    def execute(self):
        """
        Runs the requested Orac run and returns an updated
        HPC_Drug.structures.protein.Protein instance

        (The only pubblic method of this class)
        """

        #writes self.orac_in_file
        self._write_template_on_file()

        #run Orac
        self._run()

        #Orac removes the chain id from the pdb file
        #I put it back
        self._write_chain_in_pdb()

        self.Protein.pdb_file = self.output_pdb_file
 
        return self.Protein
