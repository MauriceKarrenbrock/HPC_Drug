######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains functions and classes for the creation of a gromacs input
"""
import Bio.PDB
import subprocess
import os
import shutil
import math

from HPC_Drug.structures import ligand
from HPC_Drug.structures import protein
from HPC_Drug import file_manipulation
from HPC_Drug import important_lists
from HPC_Drug import funcs4orac
from HPC_Drug import pipeline_functions
from HPC_Drug import orient
from HPC_Drug import funcs4pbs
from HPC_Drug import funcs4slurm
from HPC_Drug.auxiliary_functions import get_iterable

def residue_substitution(Protein, substitution = 'standard', ph = 7.0):
    """Takes a protein instance, and returns one

    renames the resnames of the metal binding residues in order
    to generate the right force field
    
    decides how histidines should be protonated (ph > 5 HSD, ph < 5 HIS except for the metal binding ones)"""

    #You can choose any implemented substitution
    #or create your own
    if substitution == 'standard':
        def substitute(input_list) :
            return standard_gromacs_substitutions(input_list)
    
    elif substitution == 'custom_zinc':
        def substitute(input_list) :
            return custom_zinc_gromacs_substitutions(input_list)
    
    else:
        raise NotImplementedError(substitution)

    p = Bio.PDB.PDBParser()
    struct = p.get_structure(Protein.protein_id, Protein.pdb_file)

    #I iterate through the structure
    for model in struct:
        for chain in model:
            for residue in chain:

                #search for the right residues
                res_id = str(residue.id[1])

                #check if they are bounding a metal and are not a disulfide bond
                if res_id in Protein.substitutions_dict.keys():
                    if Protein.substitutions_dict[res_id][2] in important_lists.metals:

                        #make substitutions with the selected function
                        residue.resname = substitute(Protein.substitutions_dict[res_id])

                    #if a cysteine is in a disulfide bond it is called CYX
                    elif residue.resname in important_lists.cyst_resnames:
                        residue.resname = _disulf_cysteine_gromacs_substitutions(
                                                    input_list = Protein.substitutions_dict[res_id],
                                                    resname = residue.resname
                                                    )

                #give the right name to histidines in order to get the right protonation
                elif residue.resname in important_lists.hist_resnames:
                    
                    if ph < 6.0:
                        residue.resname = 'HIS'
                    
                    else:
                        residue.resname = 'HID'


    Protein.structure = struct
    Protein.write(file_name = Protein.pdb_file, struct_type = 'biopython')

    return Protein


def standard_gromacs_substitutions(input_list):
    """It is called by residue_substitution if run with standard
    substitution option
    takes a tuple long 3 with (resname, binding atom, metal name)
    
    returns the right resname"""
    
    if input_list[0] in important_lists.hist_resnames:

        if input_list[1] == 'NE2':
            resname = 'HID'
        
        elif input_list[1] == 'ND1':
            resname = 'HIE'
        
        else:
            raise NotImplementedError(f"{input_list[1]} is not an implemented atom name for histidine")

    elif input_list[0] in important_lists.cyst_resnames:
        resname = 'CYM'

    else:
        print(input_list[0], 'substitution not implemented, going on like nothing happened')
        resname = input_list[0]

    return resname

def custom_zinc_gromacs_substitutions(input_list):
    """It is a custom version of standard_gromacs_substitutions(list) function
    with custom names for zinc binding residues"""
    
    resname = standard_gromacs_substitutions(input_list)

    if resname == 'HID' and input_list[2] == 'ZN':
        resname = 'HDZ'
    
    elif resname == 'HIE' and input_list[2] == 'ZN':
        resname = 'HEZ'

    elif resname == 'CIM' and input_list[2] == 'ZN':
        resname = 'CIZ'

    return resname

def _disulf_cysteine_gromacs_substitutions(input_list, resname = None):
    """Gives the CYX resname to the disulfide bond cysteines"""
    
    if input_list[2].lower().strip() == "disulf":

        return 'CYX'

    elif resname != None:

        return resname

    else:

        return input_list[0]


def gro2pdb(gro_file, pdb_file = None, chain = 'A', gromacs_path = 'gmx'):
    """converts a gro file in a pdb using gmx editconf and adds
    the chain identifier to the pdb (it is missing and it confuses many pdb parsers)

    gro_file :: string, the gro file to convert
    pdb_file :: string, optional, the name of the output file
    chain :: string, optional, default = A, it is the chain value that will be inserted in the pdb
    gromacs_path :: string, the gromacs executable default gmx

    returns a string of the pdb file name
    """

    def write_chain_in_pdb(chain_id = None, pdb_file = None):
        """This is a patch because orac removes the chain id from
        pdb files and this confuses some pdb parsers"""

        with open(pdb_file, 'r') as f:
            lines = f.readlines()

        with open(pdb_file, 'w') as w:

            for i in range(len(lines)):

                if lines[i][0:4] == 'ATOM' or lines[i][0:6] == 'HETATM' or lines[i][0:3] == 'TER':

                    lines[i] = lines[i][:21] + chain_id + lines[i][22:]

                w.write(f"{lines[i].strip()}\n")



    if pdb_file == None:
        pdb_file = gro_file.rsplit('.', 1)[0].strip() + ".pdb"

    subprocess.run(f"{gromacs_path} editconf -f {gro_file} -o {pdb_file}",
                        shell = True,
                        check = True)

    write_chain_in_pdb(chain_id = chain, pdb_file = pdb_file)

    return pdb_file


class GromacsInput(object):
    """Generic Gromacs input template object"""

    def __init__(self,
                input_filename = None,
                output_filename = None,
                Protein = None,
                Ligand = None,
                protein_tpg_file = None,
                solvent_model = None,
                MD_program_path = 'gmx'):
        
        self.Protein = Protein
        self.Ligand = Ligand

        #input filename is the cleaned gro file with both protein and ligand
        #but no solvent or trash HETATMS
        self.input_filename = input_filename
        
        self.output_filename = output_filename

        self.protein_tpg_file = protein_tpg_file

        self.solvent_model = solvent_model
        
        self.MD_program_path = MD_program_path
        if self.MD_program_path == None:
            raise Exception('Need a MD_program_path (for example ~/bin/gmx)')

        self.template = None

        self.command_string = None
    

    def write_template_on_file(self, template = None, filename = None):
        """Writes the objects template on {filename} file"""
        if template == None:
            template = self.template

        if filename == None:
            filename  = self.output_filename

        with open(filename, 'w') as f:
            for line in template:
                f.write(f"{line}\n")

        return filename
    

    def interact_with_gromacs(self, string = None):
        """Interacts with gromacs running string
        with subprocess.run
        string must contain the gromacs path"""

        print("Running Gromacs")
        r = subprocess.run(string,
                        shell = True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        universal_newlines=True)

        print(r.stdout)
        print(r.stderr)

        if r.returncode != 0:
            raise Exception(f"Gromacs failure\n{r.stdout}\n{r.stderr}")


    def execute(self,
                template = None,
                filename = None,
                output_file = None,
                command_string = None):

        if template == None:
            template = self.template

        if filename == None:
            filename = self.output_filename
        
        if output_file == None:
            output_file = self.input_filename
        
        if command_string == None:
            command_string = self.command_string

        filename = self.write_template_on_file(template, filename)

        self.interact_with_gromacs(string = command_string)

        try:
            #make a pdb file from the gro file
            self.Protein.pdb_file = gro2pdb(gro_file = output_file,
                                            pdb_file = self.Protein.pdb_file,
                                            chain = self.Protein.chain,
                                            gromacs_path = self.MD_program_path)

        except:
            pass
 
        return output_file


class GromacsMakeProteinGroTop(GromacsInput):
    """Takes a protein instance and creates it's .top and .gro file"""

    def __init__(self,
                input_filename = None,
                output_filename = "choices.txt",
                Protein = None,
                Ligand = None,
                protein_tpg_file = '6',
                solvent_model = '7',
                MD_program_path = None):

        super().__init__(self,
                #input_filename = input_filename,
                output_filename = output_filename,
                Protein = Protein,
                Ligand = Ligand,
                protein_tpg_file = protein_tpg_file,
                solvent_model = solvent_model)
        #protein_tpg_file is actually the protein model chosen through choices_file
        
        self.command_string = [f"{self.MD_program_path} pdb2gmx -f {self.Protein.pdb_file} -o {self.Protein.protein_id}.gro -p {self.Protein.protein_id}.top < {self.output_filename}",
            f"{self.MD_program_path} pdb2gmx -ignh -f {self.Protein.pdb_file} -o {self.Protein.protein_id}.gro -p {self.Protein.protein_id}.top < {self.output_filename}"]

    
        self.template = [
                        f"{self.protein_tpg_file}",
                        f"{self.solvent_model}"]
    
    def execute(self,
            template = None,
            filename = None,
            output_file = None,
            Protein = None):

        """Takes a protein instance and returns one"""

        if Protein == None:
            Protein = self.Protein

        #sometimes gromacs complains about some existing hydrogens, in case I tell it to ignore them -ignh
        try:
            choices_file = super().execute(template = template,
                            filename = filename,
                            output_file = output_file,
                            command_string = self.command_string[0])

        except:
            choices_file = super().execute(template = template,
                            filename = filename,
                            output_file = output_file,
                            command_string = self.command_string[1])


        choices_file = None

        Protein.gro_file = f"{Protein.protein_id}.gro"
        Protein.top_file = f"{Protein.protein_id}.top"

        return Protein




class GromacsMakeJoinedProteinLigandTopGro(GromacsInput):
    """Given the separated ligand and protein pdb (through theire instances)
    gives a joined pdb gro and top file
    returns a Protein instance containing them"""

    def execute(self,
                Protein = None,
                Ligand = None,
                MD_program_path = None):
        
        if Protein == None:
            Protein = self.Protein

        if Ligand == None:
            Ligand = self.Ligand

        if self.Ligand == None:
            return Protein
        
        if MD_program_path == None:
            MD_program_path = self.MD_program_path

        if self.output_filename == None:
            self.output_filename = f"{Protein.protein_id}_joined.pdb"
        
        #make pdb for the protein from the gro file (in order to get the hidrogens)
        self.interact_with_gromacs(string = f'{MD_program_path} editconf -f {Protein.gro_file} -o {Protein.pdb_file}')

        #join the protein and ligand pdb
        Protein = funcs4orac.join_ligand_and_protein_pdb(Protein = self.Protein,
                                                        Ligand = self.Ligand,
                                                        output_filename = self.output_filename)

        #Create the joined gro file
        self.interact_with_gromacs(string = f'{MD_program_path} editconf -f {Protein.pdb_file} -d 1 -bt triclinic -angles 90 90 90 -o {Protein.protein_id}_joined.gro')
                                                        
        Protein.gro_file = f"{Protein.protein_id}_joined.gro"

        Protein.pdb_file = self.output_filename

        #edit and rename the top file
        Protein = self._edit_top_file(Protein = Protein, Ligand = Ligand)
        
        return Protein

    def _edit_top_file(self, Protein, Ligand):
        """Adds the needed #include and other informations
        to the protein top file in order
        to include the ligand"""

        with open(Protein.top_file, 'r') as f:
            top = f.readlines()

        itp_insertion_string = ''
        for lgand in get_iterable.get_iterable(Ligand):
            itp_insertion_string = itp_insertion_string + f'#include "{lgand.itp_file}"\n'
            compound_string = f'{lgand.resname}              1'
            top.append(compound_string)

        for i in range(len(top)):
            if '[ moleculetype ]' in top[i] or '; Include chain topologies' in top[i]:
                top[i] = itp_insertion_string + '\n' + top[i]
                break
            
        with open(f'{Protein.protein_id}_joined.top', 'w') as f:
            
            for line in top:
                f.write(f"{line}\n")

        Protein.top_file = f'{Protein.protein_id}_joined.top'
        
        return Protein

class GromacsFirstOptimization(GromacsInput):

    """Makes the first optimization of the joined structure"""

    def __init__(self,
                input_filename = None,
                output_filename = None,
                Protein = None,
                Ligand = None,
                protein_tpg_file = None,
                solvent_model = None,
                MD_program_path = 'gmx'):

        super().__init__(input_filename = input_filename,
                        output_filename = output_filename,
                        Protein = Protein,
                        Ligand = Ligand,
                        protein_tpg_file = protein_tpg_file,
                        solvent_model = solvent_model,
                        MD_program_path = MD_program_path)

        if self.input_filename == None:
            self.input_filename = f'{self.Protein.protein_id}_joined_optimized.gro'

        if self.output_filename == None:
            self.output_filename = f"{self.Protein.protein_id}_first_opt.mdp"

        #creates the .tpr and then optimizes the structure
        self.command_string = f"{self.MD_program_path} grompp -f {self.output_filename} -c {Protein.gro_file} -p {Protein.top_file} -o {Protein.protein_id}_joined_optimized.tpr -maxwarn 100 && gmx mdrun -s {Protein.protein_id}_joined_optimized.tpr  -c {self.input_filename}"

        self.template = ["integrator	= steep",
                        "nsteps		= 1000",
                        "emtol		= 100",
                        "emstep		= 0.01",
                        "nstxout 	= 1",
                        "nstenergy	= 1",
                        "rlist		= 1",
                        "coulombtype	= pme",
                        "vdw-type	= cut-off",
                        "rvdw		= 1.0",
                        "rcoulomb	= 1.0",
                        "constraints	= none"]

        

class GromacsSolvBoxInput(GromacsInput): 
    """Creates and optimizes a water box around the protein"""

    def __init__(self,
                input_filename = None,
                output_filename = None,
                Protein = None,
                Ligand = None,
                protein_tpg_file = None,
                solvent_model = "amber99sb-ildn.ff/spce.itp",
                MD_program_path = 'gmx',
                box_borders = '0.8'):

        super().__init__(input_filename = input_filename,
                        output_filename = output_filename,
                        Protein = Protein,
                        Ligand = Ligand,
                        protein_tpg_file = protein_tpg_file,
                        solvent_model = solvent_model,
                        MD_program_path = MD_program_path)

        if self.output_filename == None:
            self.output_filename = f'{self.Protein.protein_id}_solv_box.mdp'

        if self.input_filename == None:
            self.input_filename = f"{self.Protein.protein_id}_solv_box.gro"
        
        self.box_borders = box_borders

        #creates the .tpr and then optimizes the structure it is actually a list of strings
        self.command_string = [
            f"{self.MD_program_path} editconf -f {self.Protein.gro_file} -d {self.box_borders} -bt triclinic -angles 90 90 90 -o {self.input_filename}",
            f"{self.MD_program_path} solvate -cp {self.input_filename} -p {self.Protein.top_file} -o {self.input_filename}",
            f"{self.MD_program_path} grompp -f {self.output_filename} -c {self.input_filename} -p {self.Protein.top_file} -maxwarn 100 -o {self.Protein.protein_id}_solv_box.tpr",
            f"{self.MD_program_path} mdrun -s {self.Protein.protein_id}_solv_box.tpr  -c {self.input_filename}"
                       ]

        self.template = [
                        "; VARIOUS PREPROCESSING OPTIONS",
                        "; Preprocessor information: use cpp syntax.",
                        "; e.g.: -I/home/joe/doe -I/home/mary/roe",
                        "include                  =",
                        "; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)",
                        "define                   =",
                        "",
                        "; RUN CONTROL PARAMETERS",
                        "integrator               = md",
                        "; Start time and timestep in ps",
                        "tinit                    = 0",
                        "dt                       = 0.0001",
                        "nsteps                   = 100000",
                        "; For exact run continuation or redoing part of a run",
                        "init-step                = 0",
                        "; Part index is updated automatically on checkpointing (keeps files separate)",
                        "simulation-part          = 1",
                        "; mode for center of mass motion removal",
                        "comm-mode                = Linear",
                        "; number of steps for center of mass motion removal",
                        "nstcomm                  = 100",
                        "; group(s) for center of mass motion removal",
                        "comm-grps                =",
                        "",
                        "; TEST PARTICLE INSERTION OPTIONS",
                        "rtpi                     = 0.05",
                        "",
                        "; OUTPUT CONTROL OPTIONS",
                        "; Output frequency for coords (x), velocities (v) and forces (f)",
                        "nstxout                  = 1000",
                        "nstvout                  = 1000",
                        "nstfout                  = 1000",
                        "; Output frequency for energies to log file and energy file",
                        "nstlog                   = 1000",
                        "nstcalcenergy            = 100",
                        "nstenergy                = 1000",
                        "; Output frequency and precision for .xtc file",
                        "nstxtcout                = 100",
                        "xtc-precision            = 1000",
                        "; This selects the subset of atoms for the .xtc file. You can",
                        "; select multiple groups. By default all atoms will be written.",
                        "xtc-grps                 =",
                        "; Selection of energy groups",
                        "energygrps               = System",
                        "",
                        "; NEIGHBORSEARCHING PARAMETERS",
                        "; cut-off scheme (group: using charge groups, Verlet: particle based cut-offs)",
                        "; nblist update frequency",
                        "cutoff-scheme            = Verlet",
                        "nstlist                  = 20",
                        "verlet-buffer-tolerance  = 0.0001",
                        "; ns algorithm (simple or grid)",
                        "ns_type                  = grid",
                        "; Periodic boundary conditions: xyz, no, xy",
                        "pbc                      = xyz",
                        "periodic-molecules       = no",
                        "; Allowed energy drift due to the Verlet buffer in kJ/mol/ps per atom,",
                        "; a value of -1 means: use rlist",
                        "; nblist cut-off",
                        "rlist                    = 1",
                        "; long-range cut-off for switched potentials",
                        "rlistlong                = -1",
                        "",
                        "; OPTIONS FOR ELECTROSTATICS AND VDW",
                        "; Method for doing electrostatics",
                        "coulombtype              = PME",
                        "rcoulomb-switch          = 0",
                        "rcoulomb                 = 1.0",
                        "; Relative dielectric constant for the medium and the reaction field",
                        "epsilon-r                = 1",
                        "epsilon-rf               = 0",
                        "; Method for doing Van der Waals",
                        "vdw-type                 = Cut-off",
                        "; cut-off lengths",
                        "rvdw-switch              = 0",
                        "rvdw                     = 1.0",
                        "; Apply long range dispersion corrections for Energy and Pressure",
                        "DispCorr                 = EnerPres",
                        "; Extension of the potential lookup tables beyond the cut-off",
                        "table-extension          = 1",
                        "; Separate tables between energy group pairs",
                        "energygrp-table          =",
                        "; Spacing for the PME/PPPM FFT grid",
                        "fourierspacing           = 0.12",
                        "; FFT grid size, when a value is 0 fourierspacing will be used",
                        "fourier-nx               = 0",
                        "fourier-ny               = 0",
                        "fourier-nz               = 0",
                        "; EWALD/PME/PPPM parameters",
                        "pme-order                = 4",
                        "ewald-rtol               = 1e-06",
                        "ewald-geometry           = 3d",
                        "epsilon-surface          =",
                        "optimize-fft             = no",
                        "",
                        "; IMPLICIT SOLVENT ALGORITHM",
                        "implicit-solvent         = No",
                        "",
                        "; OPTIONS FOR WEAK COUPLING ALGORITHMS",
                        "; Temperature coupling",
                        "tcoupl                   = v-rescale",
                        "nsttcouple               = -1",
                        "nh-chain-length          = 1",
                        "; Groups to couple separately",
                        "tc-grps                  = System",
                        "; Time constant (ps) and reference temperature (K)",
                        "tau-t                    = 0.2",
                        "ref-t                    = 298.15",
                        "; pressure coupling",
                        "pcoupl                   = Berendsen",
                        "pcoupltype               = Isotropic",
                        "nstpcouple               = -1",
                        "; Time constant (ps), compressibility (1/bar) and reference P (bar)",
                        "tau-p                    = 0.5",
                        "compressibility          = 4.6e-5",
                        "ref-p                    = 1",
                        "; Scaling of reference coordinates, No, All or COM",
                        "refcoord-scaling         = COM",
                        "",
                        "; GENERATE VELOCITIES FOR STARTUP RUN",
                        "gen-vel                  = no",
                        "gen-temp                 = 500",
                        "gen-seed                 = 173529",
                        "",
                        "; OPTIONS FOR BONDS",
                        "constraints              = none",
                        "; Type of constraint algorithm",
                        "constraint-algorithm     = Lincs",
                        "; Do not constrain the start configuration",
                        "continuation             = no",
                        "; Use successive overrelaxation to reduce the number of shake iterations",
                        "Shake-SOR                = no",
                        "; Relative tolerance of shake",
                        "shake-tol                = 0.00001",
                        "; Highest order in the expansion of the constraint coupling matrix",
                        "lincs-order              = 5",
                        "; Number of iterations in the final step of LINCS. 1 is fine for",
                        "; normal simulations, but use 2 to conserve energy in NVE runs.",
                        "; For energy minimization with constraints it should be 4 to 8.",
                        "lincs-iter               = 2",
                        "; Lincs will write a warning to the stderr if in one step a bond",
                        "; rotates over more degrees than",
                        "lincs-warnangle          = 30",
                        "; Convert harmonic bonds to morse potentials",
                        "morse                    = no"
                        ]

    def _edit_top_file(self, Protein, water_itp = None):
        """Adds the needed #include for water
        to the protein top file"""

        if water_itp == None:
            water_itp = "amber99sb-ildn.ff/spce.itp"

        with open(Protein.top_file, 'r') as f:
            top = f.readlines()

        itp_insertion_string = f'#include "{water_itp}"'

        for i in range(len(top)):
            if '[ moleculetype ]' in top[i] or '; Include chain topologies' in top[i]:
                top[i] = itp_insertion_string + '\n' + top[i]
                break


        with open(Protein.top_file, 'w') as f:
            for line in top:
                f.write(f"{line}\n")

        
        return Protein

    def execute(self,
                template = None,
                filename = None,
                output_file = None,
                command_string = None):

        if template == None:
            template = self.template

        if filename == None:
            filename = self.output_filename
        
        if output_file == None:
            output_file = self.input_filename
        
        if command_string == None:
            command_string = self.command_string

        filename = self.write_template_on_file(template, filename)

        self.Protein = self._edit_top_file(Protein = self.Protein, water_itp = self.solvent_model)
        for string in command_string:
            self.interact_with_gromacs(string = string)

        #make a pdb file from the gro file
        self.Protein.pdb_file = gro2pdb(gro_file = output_file,
                                        pdb_file = self.Protein.pdb_file,
                                        chain = self.Protein.chain,
                                        gromacs_path = self.MD_program_path)
 
        return output_file


class GromacsREMInput(GromacsInput):
    
    def __init__(self,
                input_filename = None,
                output_filename = None,
                Protein = None,
                Ligand = None,
                protein_tpg_file = None,
                solvent_model = "amber99sb-ildn.ff/spce.itp",
                MD_program_path = 'gmx',
                kind_of_processor = 'skylake',
                number_of_cores_per_node = 64,
                use_gpu = 'auto'):

        super().__init__(input_filename = input_filename,
                        output_filename = output_filename,
                        Protein = Protein,
                        Ligand = Ligand,
                        protein_tpg_file = protein_tpg_file,
                        solvent_model = solvent_model,
                        MD_program_path = MD_program_path)


        #the tpr file
        if self.input_filename == None:
            self.input_filename = f"{self.Protein.protein_id}_REM.tpr"
        
        elif self.input_filename.rsplit('.', 1)[-1].strip() != 'tpr':
            self.input_filename = self.input_filename + '.tpr'

        if self.output_filename == None:
            self.output_filename = f"{self.Protein.protein_id}_REM.mdp"

        self.kind_of_processor = kind_of_processor
        self.number_of_cores_per_node = number_of_cores_per_node

        #an instance of orient.Orient class
        self.orient = orient.Orient(self.Protein, self.Ligand)

        #gromacs has various options to use gpu
        #auto (default) that will use all the available ones automaticly
        #cpu uses no GPU even if available
        #gpu forces the use of GPU (but in case you want to use a gpu auto would be safer and more robust)
        self.use_gpu = use_gpu.lower().strip()
        if self.use_gpu not in ('auto', 'cpu', 'gpu'):
            raise ValueError(f"{self.use_gpu} is not a valid gpu option, valid options are auto cpu gpu")



        self.template = [
            "; VARIOUS PREPROCESSING OPTIONS",
            "; Preprocessor information: use cpp syntax.",
            "; e.g.: -I/home/joe/doe -I/home/mary/roe",
            "include                  =",
            "; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)",
            "define                   =",
            "",
            "; RUN CONTROL PARAMETERS",
            "integrator               = md",
            "; Start time and timestep in ps",

            self.write_TIME_TIMESTEP_string(),

            "; For exact run continuation or redoing part of a run",
            "init-step                = 0",
            "; Part index is updated automatically on checkpointing (keeps files separate)",
            "simulation-part          = 1",
            "; mode for center of mass motion removal",
            "comm-mode                = Linear",
            "; number of steps for center of mass motion removal",
            "nstcomm                  = 100",
            "; group(s) for center of mass motion removal",
            "comm-grps                =",
            "",
            "; TEST PARTICLE INSERTION OPTIONS",
            "rtpi                     = 0.05",
            "",
            "; OUTPUT CONTROL OPTIONS",
            "; Output frequency for coords (x), velocities (v) and forces (f)",
            "nstxout                  = 100000",
            "nstvout                  = 100000",
            "nstfout                  = 100000",
            "; Output frequency for energies to log file and energy file",
            "nstlog                   = 100000",
            "nstcalcenergy            = 100",
            "nstenergy                = 100000",
            "; Output frequency and precision for .xtc file",
            "nstxtcout                = 80000",
            "xtc-precision            = 1000",
            "; This selects the subset of atoms for the .xtc file. You can",
            "; select multiple groups. By default all atoms will be written.",
            "xtc-grps                 =",
            "; Selection of energy groups",
            "energygrps               = System",
            "",
            "; NEIGHBORSEARCHING PARAMETERS",
            "; cut-off scheme (group: using charge groups, Verlet: particle based cut-offs)",
            "; nblist update frequency",
            "cutoff-scheme            = Verlet",
            "nstlist                  = 20",
            "verlet-buffer-tolerance  = 0.0001",
            "; ns algorithm (simple or grid)",
            "ns_type                  = grid",
            "; Periodic boundary conditions: xyz, no, xy",
            "pbc                      = xyz",
            "periodic-molecules       = no",
            "; Allowed energy drift due to the Verlet buffer in kJ/mol/ps per atom,",
            "; a value of -1 means: use rlist",
            "; nblist cut-off",
            "rlist                    = 1.0",
            "; long-range cut-off for switched potentials",
            "rlistlong                = -1",
            "",
            "; OPTIONS FOR ELECTROSTATICS AND VDW",
            "; Method for doing electrostatics",
            "coulombtype              = PME",
            "rcoulomb-switch          = 0",
            "rcoulomb                 = 1.0",
            "; Relative dielectric constant for the medium and the reaction field",
            "epsilon-r                = 1",
            "epsilon-rf               = 0",
            "; Method for doing Van der Waals",
            "vdw-type                 = Cut-off",
            "; cut-off lengths",
            "rvdw-switch              = 0",
            "rvdw                     = 1.0",
            "; Apply long range dispersion corrections for Energy and Pressure",
            "DispCorr                 = EnerPres",
            "; Extension of the potential lookup tables beyond the cut-off",
            "table-extension          = 1",
            "; Separate tables between energy group pairs",
            "energygrp-table          =",
            "; Spacing for the PME/PPPM FFT grid",
            "fourierspacing           = 0.12",
            "; FFT grid size, when a value is 0 fourierspacing will be used",
            "fourier-nx               = 0",
            "fourier-ny               = 0",
            "fourier-nz               = 0",
            "; EWALD/PME/PPPM parameters",
            "pme-order                = 4",
            "ewald-rtol               = 1e-06",
            "ewald-geometry           = 3d",
            "epsilon-surface          =",
            "optimize-fft             = no",
            "",
            "; IMPLICIT SOLVENT ALGORITHM",
            "implicit-solvent         = No",
            "",
            "; OPTIONS FOR WEAK COUPLING ALGORITHMS",
            "; Temperature coupling",
            "tcoupl                   = v-rescale",
            "nsttcouple               = -1",
            "nh-chain-length          = 1",
            "; Groups to couple separately",
            "tc-grps                  = System",
            "; Time constant (ps) and reference temperature (K)",
            "tau-t                    = 0.2",
            "ref-t                    = 298.15",
            "; pressure coupling",
            "pcoupl                   = Berendsen",
            "pcoupltype               = Isotropic",
            "nstpcouple               = -1",
            "; Time constant (ps), compressibility (1/bar) and reference P (bar)",
            "tau-p                    = 0.5",
            "compressibility          = 4.6e-5",
            "ref-p                    = 1",
            "; Scaling of reference coordinates, No, All or COM",
            "refcoord-scaling         = COM",
            "",
            "; GENERATE VELOCITIES FOR STARTUP RUN",
            "gen-vel                  = no",
            "gen-temp                 = 500",
            "gen-seed                 = 173529",
            "",
            "; OPTIONS FOR BONDS",
            "constraints              = all-bonds",
            "; Type of constraint algorithm",
            "constraint-algorithm     = Lincs",
            "; Do not constrain the start configuration",
            "continuation             = no",
            "; Use successive overrelaxation to reduce the number of shake iterations",
            "Shake-SOR                = no",
            "; Relative tolerance of shake",
            "shake-tol                = 0.00001",
            "; Highest order in the expansion of the constraint coupling matrix",
            "lincs-order              = 5",
            "; Number of iterations in the final step of LINCS. 1 is fine for",
            "; normal simulations, but use 2 to conserve energy in NVE runs.",
            "; For energy minimization with constraints it should be 4 to 8.",
            "lincs-iter               = 2",
            "; Lincs will write a warning to the stderr if in one step a bond",
            "; rotates over more degrees than",
            "lincs-warnangle          = 30",
            "; Convert harmonic bonds to morse potentials",
            "morse                    = no"
        ]


    def get_ns_per_day(self,  Protein = None, Ligand = None, kind_of_processor = None):
        """Get's the number of ns per day that a kind_of_processor
        processor can process on the sistem"""

        if Ligand == None:
            Ligand = self.Ligand

        if Protein == None:
            Protein = self.Protein

        if kind_of_processor == None:
            kind_of_processor = self.kind_of_processor

        number_of_atoms = self.orient.get_first_last_atom_strucure(Protein = Protein, Ligand = Ligand)
        number_of_atoms = number_of_atoms[0]

        if self.use_gpu == 'cpu':
            ns_per_day = ( 15000. / number_of_atoms ) * important_lists.processor_kind_ns_per_day_15000_atoms_for_cpu_only_runs[kind_of_processor]
        
        elif self.use_gpu in ('auto', 'gpu'):
            ns_per_day = ( 15000. / number_of_atoms ) * important_lists.ns_per_day_15000_on_gpu_accellerated_architectures

        else:
            raise ValueError(f"{self.use_gpu} is not a valid gpu option, valid options are auto cpu gpu")

        return ns_per_day

    def write_TIME_TIMESTEP_string(self):
        """writes the tinit timestep and number of steps string"""

        tinit = 0
        timestep = 0.00150 #ps
        
        #ps calculated in 24 hours by one mpi_run (BATTERIES)
        number_of_steps = self.get_ns_per_day() * 1.E+3
        #number of steps (integer)
        number_of_steps = math.ceil( number_of_steps / timestep )

        string = f"tinit                    = {tinit} \ndt                       = {timestep} \nnsteps                   = {number_of_steps}\n"

        return string

    def get_BATTERIES(self, Protein = None, Ligand = None, kind_of_processor = None):
        """Get's the number of batteries for REM"""

        if Ligand == None:
            Ligand = self.Ligand

        if Protein == None:
            Protein = self.Protein

        if kind_of_processor == None:
            kind_of_processor = self.kind_of_processor

        ns_per_day = self.get_ns_per_day(Protein = Protein, Ligand = Ligand, kind_of_processor = kind_of_processor)

        BATTERIES = math.ceil( 32. / ns_per_day )

        return BATTERIES

    def _edit_top_file(self):
        """PRIVATE"""

        hot_residues = self.orient.get_hot_residues_for_rem(Protein = self.Protein, Ligand = self.Ligand, cutoff = 4.5, residue_dist = 10.0)
        hot_ids = []
        for residue in hot_residues:
            hot_ids.append(str(residue[1]).strip())

        #add the ligand resnumber
        Sub_Parser = file_manipulation.SubstitutionParser()

        for lgand in self.Ligand:
            hot_ids.append(str(Sub_Parser.get_ligand_resnum(Protein = self.Protein,
                                                    ligand_resnames = lgand.resname,
                                                    chain_model_selection = True)[0][1]))

        with open(self.Protein.top_file, 'r') as f:
            lines = f.readlines()

        #auxiliary bool variable
        is_hot_residue = False
        with open(self.Protein.top_file, 'w') as f:
            for i in range(len(lines)):

                if lines[i][0:9] == "; residue" and lines[i][9:13].strip() in hot_ids :
                    if lines[i][14:17].strip() != 'SOL' :

                        is_hot_residue = True
                        
                        f.write(f"{lines[i].rstrip()}\n")

                        continue

                elif lines[i][0:9] == "; residue" and lines[i][9:13].strip() not in hot_ids:

                    is_hot_residue = False

                    f.write(f"{lines[i].rstrip()}\n")

                    continue

                #This is a check for the last residue in the chain (could probably be done better)
                elif lines[i][0].strip() == "[":

                    is_hot_residue = False

                    f.write(f"{lines[i].rstrip()}\n")

                    continue

                if is_hot_residue and len(lines[i]) >= 17:
                    if lines[i][16] != " ":

                        lines[i] = lines[i][:17] + "_" + lines[i][18:]

                f.write(f"{lines[i].rstrip()}\n")

    def interact_with_plumed(self, string = None):
        """Interacts with plumed running string
        with subprocess.run
        string must contain the plumed path"""

        print("Running Plumed")
        r = subprocess.run(string,
                        shell = True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        universal_newlines=True)

        print(r.stdout)
        print(r.stderr)

        if r.returncode != 0:
            raise Exception(f"Plumed failure\n{r.stdout}\n{r.stderr}")

    def write_workloadmanager_inputs(self, input_file, mpi_runs, replicas_for_run = 8):
        """private method called by self.execute()"""

        file_list = []

        #making a nested list the firs dimension reppresents the mpirun
        #and the second the various directories that will be called
        #with each mpirun (gmx -multidir)
        multidir = [[] for x in range(mpi_runs)]

        for i in range(mpi_runs):
            for j in range(replicas_for_run):

                multidir[i].append(f"BATTERY{i}/scaled{j}")
 
        slurm = funcs4slurm.SlurmInput(MD_input_file = input_file,
                                    MD_input_directories = multidir,
                                    slurm_input_file = f'{self.output_filename.rsplit(".", 1)[0]}.slr',
                                    MD_program = 'gromacs',
                                    MD_calculation_type = 'rem',
                                    number_of_cores_per_node = self.number_of_cores_per_node,
                                    max_time = "24:00:00",
                                    ntasks = self.get_BATTERIES(Protein = self.Protein, Ligand = self.Ligand, kind_of_processor = self.kind_of_processor) * 8,
                                    cpus_per_task = 8,
                                    std_out = f'{self.output_filename.rsplit(".", 1)[0]}.out',
                                    std_err = f'{self.output_filename.rsplit(".", 1)[0]}.err',
                                    use_gpu = self.use_gpu)

        file_list.append(slurm.write())

        pbs = funcs4pbs.SlurmInput(MD_input_file = input_file,
                                    MD_input_directories = multidir,
                                    slurm_input_file = f'{self.output_filename.rsplit(".", 1)[0]}.pbs',
                                    MD_program = 'gromacs',
                                    MD_calculation_type = 'rem',
                                    number_of_cores_per_node = self.number_of_cores_per_node,
                                    max_time = "24:00:00",
                                    ntasks = self.get_BATTERIES(Protein = self.Protein, Ligand = self.Ligand, kind_of_processor = self.kind_of_processor) * 8,
                                    cpus_per_task = 8,
                                    std_out = f'{self.output_filename.rsplit(".", 1)[0]}.out',
                                    std_err = f'{self.output_filename.rsplit(".", 1)[0]}.err',
                                    use_gpu = self.use_gpu)

        file_list.append(pbs.write())

        return file_list

    def make_TPR_files_script(self, mpi_runs, replicas_for_run = 8):
        """private"""

        filename = "MAKE_TPR_FILES.sh"

        string = "#!/bin/bash\n \n##THIS SCRIPT CREATES THE TPR FILES RUN IT BEFORE THE WORKLOADMANAGER ONE\n"

        for i in range(mpi_runs):
            for j in range(replicas_for_run):
                string = string + f"gmx grompp -maxwarn 100 -o BATTERY{i}/scaled{j}/{self.input_filename} -f BATTERY{i}/{self.output_filename} -p BATTERY{i}/{self.Protein.protein_id}_scaled_{j}.top -c BATTERY{i}/{self.Protein.gro_file.rsplit('/', 1)[-1]} \n"


        with open(filename, 'w') as f:
            f.write(string)

        return filename

    def _get_hamiltonian_scaling_values(self):
        """
        Scales the hamiltonian with a geometrical progression
        scale(m) =scale^(m/(nprocs−1)) with scale = 0.2 and nprocs = 8
        0 <= m <= nprocs -1

        for more information check the orac manual http://www.chim.unifi.it/orac/orac-manual.pdf
        page 123 """

        nprocs = 8.0
        scale = 0.2

        #instantiating the list and putting the value for scale^0 = 1.0
        hamiltonian_scaling_values = [1.0]

        for m in range(1, 8):

            scale_m = scale ** ( m / (nprocs -1) )

            hamiltonian_scaling_values.append(scale_m)

        hamiltonian_scaling_values = tuple(hamiltonian_scaling_values)

        return hamiltonian_scaling_values


    def execute(self):
        """This method does not run gromacs but
        creates the input to make a REM simulation on a
        HPC cluster"""

        self.output_filename = self.write_template_on_file(self.template, self.output_filename)

        number_of_mpiruns = self.get_BATTERIES()

        #create the elaborated top file for plumed
        self.interact_with_gromacs(string = f"{self.MD_program_path} grompp -f {self.output_filename} -c {self.Protein.gro_file} -p {self.Protein.top_file} -maxwarn 100 -pp {self.Protein.top_file.rsplit('.', 1)[0]}_elaborated.top")
        self.Protein.top_file = f"{self.Protein.top_file.rsplit('.', 1)[0]}_elaborated.top"

        #finds the hot residues and edits the elaborated top file accordingly
        self._edit_top_file()

        #use plumed to make the scaled topologies

        hamiltonian_scaling_values = self._get_hamiltonian_scaling_values()

        for counter, value in enumerate(hamiltonian_scaling_values):
            self.interact_with_plumed(string = f"plumed partial_tempering {value} < {self.Protein.top_file} > {self.Protein.protein_id}_scaled_{counter}.top")

        #make an empty_plumed.dat file for plumed
        with open("empty_plumed.dat", 'w') as f:
            f.write(" ")

        i = 0
        j = 0
        for i in range(number_of_mpiruns):

            for j in range(len(hamiltonian_scaling_values)):

                #creates the REM directory that will be copied to the HPC cluster
                os.makedirs(f"{self.Protein.protein_id}_REM/BATTERY{i}/scaled{j}", exist_ok=True)

                #copy the dummy empty_plumed.dat file
                shutil.copy("empty_plumed.dat", f"{self.Protein.protein_id}_REM/BATTERY{i}/scaled{j}")

            #copy the needed files in the new directories
            j = 0
            for j in (self.Protein.gro_file, self.output_filename):
                shutil.copy(j, f"{self.Protein.protein_id}_REM/BATTERY{i}")

            #copy the scaled top files
            j = 0
            for j in range(len(hamiltonian_scaling_values)):
                shutil.copy(f"{self.Protein.protein_id}_scaled_{j}.top", f"{self.Protein.protein_id}_REM/BATTERY{i}")
        
        #make and copy the workload manager input files
        #write workload manager input for different workload managers (slurm pbs ...)
        workload_files = self.write_workloadmanager_inputs(input_file = self.input_filename, mpi_runs = number_of_mpiruns, replicas_for_run = len(hamiltonian_scaling_values))
        for wl_file in workload_files:
            shutil.copy(wl_file, f"{self.Protein.protein_id}_REM")


        #make and copy the script that will make the tpr files in loco
        TPR_file_script = self.make_TPR_files_script(mpi_runs = number_of_mpiruns, replicas_for_run = len(hamiltonian_scaling_values))
        shutil.copy(TPR_file_script, f"{self.Protein.protein_id}_REM")

        return self.output_filename

        


class GromacsNativeREMInput(GromacsInput):
    #WORK IN PROGRESS
    
    def __init__(self,
                input_filename = None,
                output_filename = None,
                Protein = None,
                Ligand = None,
                protein_tpg_file = None,
                solvent_model = "amber99sb-ildn.ff/spce.itp",
                MD_program_path = 'gmx',
                kind_of_processor = 'skylake',
                number_of_cores_per_node = 64,
                use_gpu = 'auto'):

        super().__init__(input_filename = input_filename,
                        output_filename = output_filename,
                        Protein = Protein,
                        Ligand = Ligand,
                        protein_tpg_file = protein_tpg_file,
                        solvent_model = solvent_model,
                        MD_program_path = MD_program_path)


        #the tpr file
        if self.input_filename == None:
            self.input_filename = f"{self.Protein.protein_id}_NativeREM.tpr"
        
        elif self.input_filename.rsplit('.', 1)[-1].strip() != 'tpr':
            self.input_filename = self.input_filename + '.tpr'

        if self.output_filename == None:
            self.output_filename = f"{self.Protein.protein_id}_NativeREM.mdp"

        self.kind_of_processor = kind_of_processor
        self.number_of_cores_per_node = number_of_cores_per_node

        #an instance of orient.Orient class
        self.orient = orient.Orient(self.Protein, self.Ligand)

        #gromacs has various options to use gpu
        #auto (default) that will use all the available ones automaticly
        #cpu uses no GPU even if available
        #gpu forces the use of GPU (but in case you want to use a gpu auto would be safer and more robust)
        self.use_gpu = use_gpu.lower().strip()
        if self.use_gpu not in ('auto', 'cpu', 'gpu'):
            raise ValueError(f"{self.use_gpu} is not a valid gpu option, valid options are auto cpu gpu")



        self.template = [
            "; VARIOUS PREPROCESSING OPTIONS",
            "; Preprocessor information: use cpp syntax.",
            "; e.g.: -I/home/joe/doe -I/home/mary/roe",
            "include                  =",
            "; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)",
            "define                   =",
            "",
            "; RUN CONTROL PARAMETERS",
            "integrator               = md",
            "; Start time and timestep in ps",

            self.write_TIME_TIMESTEP_string(),

            "; For exact run continuation or redoing part of a run",
            "init-step                = 0",
            "; Part index is updated automatically on checkpointing (keeps files separate)",
            "simulation-part          = 1",
            "; mode for center of mass motion removal",
            "comm-mode                = Linear",
            "; number of steps for center of mass motion removal",
            "nstcomm                  = 100",
            "; group(s) for center of mass motion removal",
            "comm-grps                =",
            "",
            "; TEST PARTICLE INSERTION OPTIONS",
            "rtpi                     = 0.05",
            "",
            "; OUTPUT CONTROL OPTIONS",
            "; Output frequency for coords (x), velocities (v) and forces (f)",
            "nstxout                  = 100000",
            "nstvout                  = 100000",
            "nstfout                  = 100000",
            "; Output frequency for energies to log file and energy file",
            "nstlog                   = 100000",
            "nstcalcenergy            = 100",
            "nstenergy                = 100000",
            "; Output frequency and precision for .xtc file",
            "nstxtcout                = 80000",
            "xtc-precision            = 1000",
            "; This selects the subset of atoms for the .xtc file. You can",
            "; select multiple groups. By default all atoms will be written.",
            "xtc-grps                 =",
            "; Selection of energy groups",
            "energygrps               = System",
            "",
            "; NEIGHBORSEARCHING PARAMETERS",
            "; cut-off scheme (group: using charge groups, Verlet: particle based cut-offs)",
            "; nblist update frequency",
            "cutoff-scheme            = Verlet",
            "nstlist                  = 20",
            "verlet-buffer-tolerance  = 0.0001",
            "; ns algorithm (simple or grid)",
            "ns_type                  = grid",
            "; Periodic boundary conditions: xyz, no, xy",
            "pbc                      = xyz",
            "periodic-molecules       = no",
            "; Allowed energy drift due to the Verlet buffer in kJ/mol/ps per atom,",
            "; a value of -1 means: use rlist",
            "; nblist cut-off",
            "rlist                    = 1.0",
            "; long-range cut-off for switched potentials",
            "rlistlong                = -1",
            "",
            "; OPTIONS FOR ELECTROSTATICS AND VDW",
            "; Method for doing electrostatics",
            "coulombtype              = PME",
            "rcoulomb-switch          = 0",
            "rcoulomb                 = 1.0",
            "; Relative dielectric constant for the medium and the reaction field",
            "epsilon-r                = 1",
            "epsilon-rf               = 0",
            "; Method for doing Van der Waals",
            "vdw-type                 = Cut-off",
            "; cut-off lengths",
            "rvdw-switch              = 0",
            "rvdw                     = 1.0",
            "; Apply long range dispersion corrections for Energy and Pressure",
            "DispCorr                 = EnerPres",
            "; Extension of the potential lookup tables beyond the cut-off",
            "table-extension          = 1",
            "; Separate tables between energy group pairs",
            "energygrp-table          =",
            "; Spacing for the PME/PPPM FFT grid",
            "fourierspacing           = 0.12",
            "; FFT grid size, when a value is 0 fourierspacing will be used",
            "fourier-nx               = 0",
            "fourier-ny               = 0",
            "fourier-nz               = 0",
            "; EWALD/PME/PPPM parameters",
            "pme-order                = 4",
            "ewald-rtol               = 1e-06",
            "ewald-geometry           = 3d",
            "epsilon-surface          =",
            "optimize-fft             = no",
            "",
            "; IMPLICIT SOLVENT ALGORITHM",
            "implicit-solvent         = No",
            "",
            "; OPTIONS FOR WEAK COUPLING ALGORITHMS",
            "; Temperature coupling",
            "tcoupl                   = v-rescale",
            "nsttcouple               = -1",
            "nh-chain-length          = 1",
            "; Groups to couple separately",
            "tc-grps                  = System",
            "; Time constant (ps) and reference temperature (K)",
            "tau-t                    = 0.2",
            "ref-t                    = 298.15",
            "; pressure coupling",
            "pcoupl                   = Berendsen",
            "pcoupltype               = Isotropic",
            "nstpcouple               = -1",
            "; Time constant (ps), compressibility (1/bar) and reference P (bar)",
            "tau-p                    = 0.5",
            "compressibility          = 4.6e-5",
            "ref-p                    = 1",
            "; Scaling of reference coordinates, No, All or COM",
            "refcoord-scaling         = COM",
            "",
            "; GENERATE VELOCITIES FOR STARTUP RUN",
            "gen-vel                  = no",
            "gen-temp                 = 500",
            "gen-seed                 = 173529",
            "",
            "; OPTIONS FOR BONDS",
            "constraints              = all-bonds",
            "; Type of constraint algorithm",
            "constraint-algorithm     = Lincs",
            "; Do not constrain the start configuration",
            "continuation             = no",
            "; Use successive overrelaxation to reduce the number of shake iterations",
            "Shake-SOR                = no",
            "; Relative tolerance of shake",
            "shake-tol                = 0.00001",
            "; Highest order in the expansion of the constraint coupling matrix",
            "lincs-order              = 5",
            "; Number of iterations in the final step of LINCS. 1 is fine for",
            "; normal simulations, but use 2 to conserve energy in NVE runs.",
            "; For energy minimization with constraints it should be 4 to 8.",
            "lincs-iter               = 2",
            "; Lincs will write a warning to the stderr if in one step a bond",
            "; rotates over more degrees than",
            "lincs-warnangle          = 30",
            "; Convert harmonic bonds to morse potentials",
            "morse                    = no"
        ]


    def get_ns_per_day(self,  Protein = None, Ligand = None, kind_of_processor = None):
        """Get's the number of ns per day that a kind_of_processor
        processor can process on the sistem"""

        if Ligand == None:
            Ligand = self.Ligand

        if Protein == None:
            Protein = self.Protein

        if kind_of_processor == None:
            kind_of_processor = self.kind_of_processor

        number_of_atoms = self.orient.get_first_last_atom_strucure(Protein = Protein, Ligand = Ligand)
        number_of_atoms = number_of_atoms[0]

        if self.use_gpu == 'cpu':
            ns_per_day = ( 15000. / number_of_atoms ) * important_lists.processor_kind_ns_per_day_15000_atoms_for_cpu_only_runs[kind_of_processor]
        
        elif self.use_gpu in ('auto', 'gpu'):
            ns_per_day = ( 15000. / number_of_atoms ) * important_lists.ns_per_day_15000_on_gpu_accellerated_architectures

        else:
            raise ValueError(f"{self.use_gpu} is not a valid gpu option, valid options are auto cpu gpu")

        return ns_per_day

    def write_TIME_TIMESTEP_string(self):
        """writes the tinit timestep and number of steps string"""

        tinit = 0
        timestep = 0.00150 #ps
        
        #ps calculated in 24 hours by one mpi_run (BATTERIES)
        number_of_steps = self.get_ns_per_day() * 1.E+3
        #number of steps (integer)
        number_of_steps = math.ceil( number_of_steps / timestep )

        string = f"tinit                    = {tinit} \ndt                       = {timestep} \nnsteps                   = {number_of_steps}\n"

        return string

    def get_BATTERIES(self, Protein = None, Ligand = None, kind_of_processor = None):
        """Get's the number of batteries for REM"""

        if Ligand == None:
            Ligand = self.Ligand

        if Protein == None:
            Protein = self.Protein

        if kind_of_processor == None:
            kind_of_processor = self.kind_of_processor

        ns_per_day = self.get_ns_per_day(Protein = Protein, Ligand = Ligand, kind_of_processor = kind_of_processor)

        BATTERIES = math.ceil( 32. / ns_per_day )

        return BATTERIES

    def _edit_top_file(self):
        """PRIVATE"""

        hot_residues = self.orient.get_hot_residues_for_rem(Protein = self.Protein, Ligand = self.Ligand, cutoff = 4.5, residue_dist = 10.0)
        hot_ids = []
        for residue in hot_residues:
            hot_ids.append(str(residue[1]).strip())

        #add the ligand resnumber
        Sub_Parser = file_manipulation.SubstitutionParser()

        for lgand in self.Ligand:
            hot_ids.append(str(Sub_Parser.get_ligand_resnum(Protein = self.Protein,
                                                    ligand_resnames = lgand.resname,
                                                    chain_model_selection = True)[0][1]))

        with open(self.Protein.top_file, 'r') as f:
            lines = f.readlines()

        #auxiliary bool variable
        is_hot_residue = False
        with open(self.Protein.top_file, 'w') as f:
            for i in range(len(lines)):

                if lines[i][0:9] == "; residue" and lines[i][9:13].strip() in hot_ids :
                    if lines[i][14:17].strip() != 'SOL' :

                        is_hot_residue = True
                        
                        f.write(f"{lines[i].rstrip()}\n")

                        continue

                elif lines[i][0:9] == "; residue" and lines[i][9:13].strip() not in hot_ids:

                    is_hot_residue = False

                    f.write(f"{lines[i].rstrip()}\n")

                    continue

                #This is a check for the last residue in the chain (could probably be done better)
                elif lines[i][0].strip() == "[":

                    is_hot_residue = False

                    f.write(f"{lines[i].rstrip()}\n")

                    continue

                if is_hot_residue and len(lines[i]) >= 17:
                    if lines[i][16] != " ":

                        lines[i] = lines[i][:17] + "_" + lines[i][18:]

                f.write(f"{lines[i].rstrip()}\n")

    def interact_with_plumed(self, string = None):
        """Interacts with plumed running string
        with subprocess.run
        string must contain the plumed path"""

        print("Running Plumed")
        r = subprocess.run(string,
                        shell = True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        universal_newlines=True)

        print(r.stdout)
        print(r.stderr)

        if r.returncode != 0:
            raise Exception(f"Plumed failure\n{r.stdout}\n{r.stderr}")

    def write_workloadmanager_inputs(self, input_file, mpi_runs, replicas_for_run = 8):
        """private method called by self.execute()"""

        file_list = []

        #making a nested list the firs dimension reppresents the mpirun
        #and the second the various directories that will be called
        #with each mpirun (gmx -multidir)
        multidir = [[] for x in range(mpi_runs)]

        for i in range(mpi_runs):
            for j in range(replicas_for_run):

                multidir[i].append(f"BATTERY{i}/scaled{j}")
 
        slurm = funcs4slurm.SlurmInput(MD_input_file = input_file,
                                    MD_input_directories = multidir,
                                    slurm_input_file = f'{self.output_filename.rsplit(".", 1)[0]}.slr',
                                    MD_program = 'gromacs',
                                    MD_calculation_type = 'rem',
                                    number_of_cores_per_node = self.number_of_cores_per_node,
                                    max_time = "24:00:00",
                                    ntasks = self.get_BATTERIES(Protein = self.Protein, Ligand = self.Ligand, kind_of_processor = self.kind_of_processor) * 8,
                                    cpus_per_task = 8,
                                    std_out = f'{self.output_filename.rsplit(".", 1)[0]}.out',
                                    std_err = f'{self.output_filename.rsplit(".", 1)[0]}.err',
                                    use_gpu = self.use_gpu,
                                    plumed = False)

        file_list.append(slurm.write())

        pbs = funcs4pbs.SlurmInput(MD_input_file = input_file,
                                    MD_input_directories = multidir,
                                    slurm_input_file = f'{self.output_filename.rsplit(".", 1)[0]}.pbs',
                                    MD_program = 'gromacs',
                                    MD_calculation_type = 'rem',
                                    number_of_cores_per_node = self.number_of_cores_per_node,
                                    max_time = "24:00:00",
                                    ntasks = self.get_BATTERIES(Protein = self.Protein, Ligand = self.Ligand, kind_of_processor = self.kind_of_processor) * 8,
                                    cpus_per_task = 8,
                                    std_out = f'{self.output_filename.rsplit(".", 1)[0]}.out',
                                    std_err = f'{self.output_filename.rsplit(".", 1)[0]}.err',
                                    use_gpu = self.use_gpu,
                                    plumed = False)

        file_list.append(pbs.write())

        return file_list

    def make_TPR_files_script(self, mpi_runs, replicas_for_run = 8):
        """private"""

        filename = "MAKE_TPR_FILES_Native.sh"

        string = "#!/bin/bash\n \n##THIS SCRIPT CREATES THE TPR FILES RUN IT BEFORE THE WORKLOADMANAGER ONE\n"

        for i in range(mpi_runs):
            for j in range(replicas_for_run):
                string = string + f"gmx grompp -maxwarn 100 -o BATTERY{i}/scaled{j}/{self.input_filename} -f BATTERY{i}/Native_REM{j}.mdp -p BATTERY{i}/{self.Protein.protein_id}_scaled_{j}.top -c BATTERY{i}/{self.Protein.gro_file.rsplit('/', 1)[-1]} \n"


        with open(filename, 'w') as f:
            f.write(string)

        return filename

    def _get_hamiltonian_scaling_values(self):
        """
        Scales the hamiltonian with a geometrical progression
        scale(m) =scale^(m/(nprocs−1)) with scale = 0.2 and nprocs = 8
        0 <= m <= nprocs -1

        for more information check the orac manual http://www.chim.unifi.it/orac/orac-manual.pdf
        page 123 """

        nprocs = 8.0
        scale = 0.2

        #instantiating the list and putting the value for scale^0 = 1.0
        hamiltonian_scaling_values = [1.0]

        for m in range(1, 8):

            scale_m = scale ** ( m / (nprocs -1) )

            hamiltonian_scaling_values.append(scale_m)

        hamiltonian_scaling_values = tuple(hamiltonian_scaling_values)

        return hamiltonian_scaling_values



    def execute(self):
        """This method does not run gromacs but
        creates the input to make a REM simulation on a
        HPC cluster"""

        for enum, i in enumerate((298.15, 298.16, 298.17, 298.18, 298.19, 298.20, 298.21, 298.22)):

            for j in range(len(self.template)):
                
                if self.template[j].strip() == "":
                    pass
                elif self.template[j].strip().split()[0] == "ref-t":
                    
                    self.template[j] = f"ref-t                    = {i}"

            with open(f"Native_REM{enum}.mdp", "w") as w:

                for line in self.template:

                    w.write(line + '\n')

        number_of_mpiruns = self.get_BATTERIES()

        #create the elaborated top file for plumed
        #self.interact_with_gromacs(string = f"{self.MD_program_path} grompp -f {self.output_filename} -c {self.Protein.gro_file} -p {self.Protein.top_file} -maxwarn 100 -pp {self.Protein.top_file.rsplit('.', 1)[0]}_elaborated.top")
        #self.Protein.top_file = f"{self.Protein.top_file.rsplit('.', 1)[0]}_elaborated.top"

        #finds the hot residues and edits the elaborated top file accordingly
        self._edit_top_file()

        #use plumed to make the scaled topologies

        hamiltonian_scaling_values = self._get_hamiltonian_scaling_values()

        for counter, value in enumerate(hamiltonian_scaling_values):
            self.interact_with_plumed(string = f"plumed partial_tempering {value} < {self.Protein.top_file} > {self.Protein.protein_id}_scaled_{counter}.top")


        i = 0
        j = 0
        for i in range(number_of_mpiruns):

            for j in range(len(hamiltonian_scaling_values)):

                #creates the REM directory that will be copied to the HPC cluster
                os.makedirs(f"{self.Protein.protein_id}_NativeREM/BATTERY{i}/scaled{j}", exist_ok=True)

                

            #copy the needed files in the new directories
            j = 0
            for j in (self.Protein.gro_file, "Native_REM0.mdp", "Native_REM1.mdp", "Native_REM2.mdp", "Native_REM3.mdp", "Native_REM4.mdp",\
                "Native_REM5.mdp", "Native_REM6.mdp", "Native_REM7.mdp"):
                shutil.copy(j, f"{self.Protein.protein_id}_NativeREM/BATTERY{i}")

            #copy the scaled top files
            j = 0
            for j in range(len(hamiltonian_scaling_values)):
                shutil.copy(f"{self.Protein.protein_id}_scaled_{j}.top", f"{self.Protein.protein_id}_NativeREM/BATTERY{i}")
        
        #make and copy the workload manager input files
        #write workload manager input for different workload managers (slurm pbs ...)
        workload_files = self.write_workloadmanager_inputs(input_file = self.input_filename, mpi_runs = number_of_mpiruns, replicas_for_run = len(hamiltonian_scaling_values))
        for wl_file in workload_files:
            shutil.copy(wl_file, f"{self.Protein.protein_id}_NativeREM")


        #make and copy the script that will make the tpr files in loco
        TPR_file_script = self.make_TPR_files_script(mpi_runs = number_of_mpiruns, replicas_for_run = len(hamiltonian_scaling_values))
        shutil.copy(TPR_file_script, f"{self.Protein.protein_id}_NativeREM")

        return self.output_filename






