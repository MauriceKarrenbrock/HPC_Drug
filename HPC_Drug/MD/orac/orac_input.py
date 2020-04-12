######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
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
from HPC_Drug import orient
from HPC_Drug import important_lists
from HPC_Drug import lib
from HPC_Drug.auxiliary_functions import get_iterable
from HPC_Drug.auxiliary_functions import path
from HPC_Drug.files_IO import write_on_files
from HPC_Drug.files_IO import read_file

import os
import importlib_resources
import subprocess

class OracInput(object):
    """Super class of any Orac input template object"""

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

        self.template = []

        #an instance of orient.Orient class
        self.orient = orient.Orient(self.Protein, self.Protein.get_ligand_list())

    def _write_box(self):
        """
        private

        Writes the string about the box size
        and rotates the strucure in a smart way
        """

        self.Protein.structure = self.orient.base_change_structure()

        self.Protein.write(file_name = self.Protein.pdb_file, struct_type = 'biopython')
        
        lx, ly, lz = self.orient.create_box(self.Protein.structure)

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
        r = self.Protein.structure.get_residues()

        for i in r:
            resnum = i.id[1]
            break
        
        cutoff = 1 - resnum
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

        lx, ly, lz = self.orient.create_box(self.Protein.structure)

        nx, ny, nz = self.orient.create_solvent_grid(lx, ly, lz)

        string = f"   GENERATE   {nx}  {ny}  {nz}  !! the grid depends on the BOX"

        return string

    def _write_EWALD_PME(self):
        """
        private
        
        Writes the informations about the grid
        in the reciprocal reticle
        """
        
        lx, ly, lz = self.orient.create_box(self.Protein.structure)

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


        #########################################################################################
        #DA RIFARE

        #get the separated structures of the protein and the ligands
        self.Protein.structure, ligand_structures = self.orient.separate_protein_ligand()

        #Make sure to save the structures of each ligand
        #may come in handy
        for lig in get_iterable.get_iterable(self.orient.Ligand):
            self.Protein.add_ligand(lig)

        ############################################################################################
        
        #Get the atom number range of each segment
        protein_atoms, ligand_atoms = self.orient.protein_ligand_atom_numbers()

        string = []

        #get the distance of the centers of mass of the protein and the ligands
        #and create the string
        for i, ligand in enumerate(get_iterable.get_iterable(ligand_structures)):

            COM_Protein, COM_ligand, distance = self.orient.center_mass_distance(self.Protein.structure, ligand)
            
            #These are useless
            COM_Protein = None
            COM_ligand = None

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

    def write_LINKED_CELL(self):
        """writes the LINKED CELL string"""

        lx, ly, lz = self.orient.create_box(self.Protein.structure)

        string = f"   LINKED_CELL   {int(lx/3.0)}  {int(ly/3.0)}  {int(lz/3.0)}"

        return string

    def write_template_on_file(self):
        """Writes the objects template on {filename} file"""

        for i in range(len(self.template)):

            self.template[i] = self.template[i] + '\n'
        
        write_on_files.write_file(lines = self.template, file_name = self.orac_in_file)


    def write_chain_in_pdb(self, chain_id = None, pdb_file = None):
        """This is a patch because orac removes the chain id from
        pdb files and this confuses some pdb parsers"""

        if chain_id == None:
            chain_id = self.Protein.chain

        chain_id = chain_id.upper().strip()

        if pdb_file == None:
            pdb_file = self.output_pdb_file

        lines = read_file.read_file(file_name = pdb_file)

        for i in range(len(lines)):

            if lines[i][0:4] == 'ATOM' or lines[i][0:6] == 'HETATM' or lines[i][0:3] == 'TER':

                lines[i] = lines[i][:20] + "{0:>2}".format(chain_id) + lines[i][22:]

                lines[i] = lines[i].strip('\n') + '\n'

        write_on_files.write_file(lines = lines, file_name = pdb_file)

    def _run(self):
        """
        private
        """

        print("Running Orac")

        #Orac reads the input from the stdin
        with open(self.orac_in_file, 'r') as out_file:
            r = subprocess.run([self.MD_program_path],
                            shell = False,
                            stdout = subprocess.PIPE,
                            stderr = subprocess.PIPE,
                            stdin = out_file,
                            universal_newlines=False)

        print(r.stdout)
        print(r.stderr)

        if r.returncode != 0:
            raise Exception(f"Orac failure\n{r.stdout}\n{r.stderr}")

    def execute(self):
        """
        Runs the requested Orac run and restuns an updated
        HPC_Drug.structures.protein.Protein instance
        """

        #writes self.orac_in_file
        self.write_template_on_file()

        #run Orac
        self._run()

        #Orac removes the chain id from the pdb file
        #I put it back
        self.write_chain_in_pdb()

        self.Protein.pdb_file = self.output_pdb_file
 
        return self.Protein
