#containing functions and classes for the creation of a gromacs input
import Bio.PDB
import subprocess

from HPC_Drug import structures
from HPC_Drug import file_manipulation
from HPC_Drug import important_lists
from HPC_Drug import funcs4orac
from HPC_Drug import pipeline_functions

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
    struct = p.get_structure(Protein.protein_id, Protein.filename)

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
    Protein.filename = Protein.write_PDB(Protein.filename, 'biopython')

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

        #input filename is the cleaned pdb with both protein and ligand
        #but no solvent or trash HETATMS
        #if not explicitly given it is considered Protein.filename (standard behaviour)
        self.input_filename = input_filename
        if self.input_filename == None:
            self.input_filename = self.Protein.filename
        
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
        
        self.command_string = f"{self.MD_program_path} pdb2gmx -f {self.Protein.filename} -o {self.Protein.protein_id}.gro -p {self.Protein.protein_id}.top < {self.output_filename}"

    
        self.template = [
                        f"{self.protein_tpg_file}",
                        f"{self.solvent_model}"]
    
    def execute(self,
            template = None,
            filename = None,
            output_file = None,
            command_string = None,
            Protein = None):

        """Takes a protein instance and returns one"""

        if Protein == None:
            Protein = self.Protein

        choices_file = super().execute(template = template,
                        filename = filename,
                        output_file = output_file,
                        command_string = command_string)
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
        
        if MD_program_path == None:
            MD_program_path = self.MD_program_path

        if self.output_filename == None:
            self.output_filename = f"{Protein.protein_id}_joined.pdb"
        
        #make pdb for the protein from the gro file (in order to get the hidrogens)
        self.interact_with_gromacs(string = f'{MD_program_path} editconf -f {Protein.gro_file} -o {Protein.filename}')

        #join the protein and ligand pdb
        Protein = funcs4orac.join_ligand_and_protein_pdb(Protein = self.Protein,
                                                        Ligand = self.Ligand,
                                                        output_filename = self.output_filename)

        #Create the joined gro file
        self.interact_with_gromacs(string = f'{MD_program_path} editconf -f {Protein.filename} -d 1 -bt triclinic -angles 90 90 90 -o {Protein.protein_id}_joined.gro')
                                                        
        Protein.gro_file = f"{Protein.protein_id}_joined.gro"

        Protein.filename = self.output_filename

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
        for ligand in pipeline_functions.get_iterable(Ligand):
            itp_insertion_string = itp_insertion_string + f'#include "{ligand.itp_file}"\n'
            compound_string = f'{ligand.ligand_resname}              1'
            top.append(compound_string)

        for i in range(len(top)):
            if '[ moleculetype ]' in top[i]:
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
                box_borders = '1'):

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
                        "nsteps                   = 5000000",
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
            if '[ moleculetype ]' in top[i]:
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
 
        return output_file
