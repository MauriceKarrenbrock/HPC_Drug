#contains all the classes and functions needed fot the orac MD program
from HPC_Drug import structures
import Bio.PDB
from HPC_Drug import file_manipulation
from HPC_Drug import important_lists
from HPC_Drug import pipeline_functions
from HPC_Drug import lib
import importlib_resources
import subprocess
import prody
from HPC_Drug import orient

def residue_substitution(Protein, substitution = 'standard', ph = 7.0):
    """Takes a protein instance, and returns one

    renames the resnames of the metal binding residues in order
    to generate the right force field
    
    decides how histidines should be protonated (ph > 5 HSD, ph < 5 HIS except for the metal binding ones)"""

    #You can choose any implemented substitution
    #or create your own
    if substitution == 'standard':
        def substitute(input_list) :
            return standard_orac_substitutions(input_list)
    
    elif substitution == 'custom_zinc':
        def substitute(input_list) :
            return custom_zinc_orac_substitutions(input_list)
    
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

                #give the right name to histidines in order to get the right protonation
                elif residue.resname in important_lists.hist_resnames:
                    
                    if ph < 6.0:
                        residue.resname = 'HIS'
                    
                    else:
                        residue.resname = 'HSD'

    Protein.structure = struct
    Protein.filename = Protein.write_PDB(Protein.filename, 'biopython')

    return Protein

def standard_orac_substitutions(input_list):
    """It is called by residue_substitution if run with standard
    substitution option
    takes a tuple long 3 with (resname, binding atom, metal name)
    
    returns the right resname"""
    
    if input_list[0] in important_lists.hist_resnames:

        if input_list[1] == 'NE2':
            resname = 'HSD'
        
        elif input_list[1] == 'ND1':
            resname = 'HSE'
        
        else:
            raise NotImplementedError(f"{input_list[1]} is not an implemented atom name for histidine")

    elif input_list[0] in important_lists.cyst_resnames:
        #This isn't the right name to write in the input seqres for Orac
        #But as Biopython strucures can only contain 3 letter resnames
        #I am giving this as a temporary name
        #It will be changed to CYSM in the get seqres function
        #This part shall be redone better
        resname = 'CYM'

    else:
        print(input_list[0], 'substitution not implemented, going on like nothing happened')
        resname = input_list[0]

    return resname

def custom_zinc_orac_substitutions(input_list):
    """It is a custom version of standard_gromacs_substitutions(list) function
    with custom names for zinc binding residues"""
    
    resname = standard_orac_substitutions(input_list)

    if resname == 'HSD' and input_list[2] == 'ZN':
        resname = 'HDZ'
    
    elif resname == 'HSE' and input_list[2] == 'ZN':
        resname = 'HEZ'

    elif resname == 'CYM' and input_list[2] == 'ZN':
        resname = 'CYZ'

    return resname

def join_ligand_and_protein_pdb(Protein = None, Ligand = None, output_filename = None):
    """Puts the ligand and the protein together again in a pdb
    
    returns a Protein istance with a updated Protein.filename"""

    if Protein == None:
        raise Exception('Protein cannot be None type')

    if Ligand == None:
        print("There is no ligand to join, will go on pretending nothing happened")
        return Protein
    
    if output_filename == None:
        output_filename = f"{Protein.protein_id}_joined.pdb"

    # A very bad way to write first the protein atoms and then the ligands ones
    #on the selected pdb file
    #Shall implement something better in the future
    with open(output_filename, 'w') as joined:
                
        with open(Protein.filename, 'r') as prot:
            for line in prot:
                if 'ATOM' in line or 'HETATM' in line:
 
                    joined.write(line)

        for ligand in pipeline_functions.get_iterable(Ligand):
            with open(ligand.ligand_filename, 'r') as lig:
                for line in lig:
                    if 'ATOM' in line or 'HETATM' in line:
   
                        joined.write(line)

    Protein.filename = output_filename

    return Protein



class OracInput(object):
    """Generic Orac input template object"""

    def __init__(self,
                input_filename = None,
                output_filename = None,
                Protein = None,
                Ligand = None,
                protein_tpg_file = None,
                protein_prm_file = None,
                solvent_pdb = None,
                MD_program_path = None):
        
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
        #if no path is given searches the standard amber tpg inside lib module
        if self.protein_tpg_file == None:
            with importlib_resources.path('HPC_Drug.lib', 'amber99sb-ildn.tpg') as path:
                self.protein_tpg_file = path#.resolve()
            
        
        self.protein_prm_file = protein_prm_file
        #if no path is given searches the standard amber prm inside lib module
        if self.protein_prm_file == None:
            with importlib_resources.path('HPC_Drug.lib', 'amber99sb-ildn.prm') as path:
                self.protein_prm_file = path#.resolve()

        self.solvent_pdb = solvent_pdb
        #if no path is given searches the standard water.pdb inside lib module
        if self.solvent_pdb == None:
            with importlib_resources.path('HPC_Drug.lib', 'water.pdb') as path:
                self.solvent_pdb = path#.resolve()
        
        self.MD_program_path = MD_program_path
        if self.MD_program_path == None:
            raise Exception('Need a MD_program_path (for example ~/ORAC/trunk/bin/orac)')

        self.template = None

        #an instance of orient.Orient class
        self.orient = orient.Orient(self.Protein, self.Ligand)

    def write_box(self):
        """Writes the string about the box size
        and rotates the strucure in a smart way"""

        self.Protein.structure = self.orient.base_change_structure()

        self.Protein.filename = self.Protein.write_PDB(filename = self.Protein.filename,
                                                    struct_type = 'biopython')
        
        lx, ly, lz = self.orient.create_box(self.Protein.structure)

        string = f"   CRYSTAL    {lx:.2f}    {ly:.2f}    {lz:.2f}  !! simulation BOX"

        return string

    def write_sulf_bond_string(self, Protein = None):
        """Writes the part of the template inherent to the disulf bonds"""

        if Protein == None:
            Protein = self.Protein

        cutoff = self._get_protein_resnumber_cutoff(Protein = Protein)
        
        tmp_sulf = []
        
        for bond in Protein.sulf_bonds:

            #Correcting with the right cutoff
            bound_1 = int(bond[0].strip()) + cutoff
            bound_2 = int(bond[0].strip()) + cutoff

            string = f"   bond 1sg 2sg residue   {bound_1}  {bound_2}"
            tmp_sulf.append(string)
        
        string = '\n'.join(tmp_sulf)

        return string

    def _get_protein_resnumber_cutoff(self, Protein = None):
        """Orac starts the residue count from one
        but pdb files often don't
        so I get the id of the resnumber of the first residue
        
        returns the cutoff to apply to start from one"""

        p = Bio.PDB.PDBParser()
        s = p.get_structure(Protein.protein_id, Protein.filename)
        r = s.get_residues()

        for i in r:
            resnum = i.id[1]
            break
        
        cutoff = 1 - resnum
        return cutoff

    def write_ligand_tpg_path(self, Ligand = None):
        """Writes the tpg file for any given ligand"""

        if Ligand == None:
            Ligand = self.Ligand
            if Ligand == None:
                return ''

        tmp = []
        for ligand in pipeline_functions.get_iterable(Ligand):

            string = f"   READ_TPG_ASCII {ligand.topology_file} !! ligand"

            tmp.append(string)

        string = '\n'.join(tmp)

        return string

    def write_ligand_prm_path(self, Ligand = None):
        """Writes the prm file for any given ligand"""

        if Ligand == None:
            Ligand = self.Ligand
            if Ligand == None:
                return ''

        tmp = []
        for ligand in pipeline_functions.get_iterable(Ligand):

            string = f"   READ_PRM_ASCII {ligand.param_file}  !! ligand"

            tmp.append(string)

        string = '\n'.join(tmp)

        return string

    def get_ligand_name_from_tpg(self, Ligand = None):
        """Gets the name of the ligand from the tpg file
        for any given ligand"""

        if Ligand == None:
            Ligand = self.Ligand
            if Ligand == None:
                return ''

        residue_strings = []

        for ligand in pipeline_functions.get_iterable(Ligand):
            with open(ligand.topology_file, 'r') as f:

                for line in f:

                    if 'RESIDUE' in line:
                        
                        string = f"      {line.split()[1].strip()} !! ligand name in tpg file"

                        residue_strings.append(string)

                        break

                else:
                    raise Exception(f'Could not find the residue name in {Ligand.topology_file}') 

        return '\n'.join(residue_strings)

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

    def execute(self, template = None, filename = None, output_pdb_file = None, MD_program_path = None):

        if MD_program_path == None:
            MD_program_path = self.MD_program_path

        if template == None:
            template = self.template

        if filename == None:
            filename = self.output_filename
        
        if output_pdb_file == None:
            output_pdb_file = self.Protein.filename.rsplit('.', 1)[0].strip() + '_optimized.pdb'

        filename = self.write_template_on_file(template, filename)

        print("Running Orac")
        r = subprocess.run([f'{MD_program_path} < {filename}'],
                        shell = True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        universal_newlines=True)

        print(r.stdout)
        print(r.stderr)

        if r.returncode != 0:
            raise Exception(f"Orac failure\n{r.stdout}\n{r.stderr}")
 
        return output_pdb_file  





class OracFirstOptimization(OracInput):
    """Takes the input file for the first optimization
    of the protein + ligand function (the 2 or more pdb files must already be mixed)
    
    Takes a Protein and any ligand instance
    the ligand needs an already calculated .tpg and .prm file

    makes the first optimization"""

    def __init__(self,
                output_filename = None,
                Protein = None,
                Ligand = None,
                protein_tpg_file = None,
                protein_prm_file = None,
                MD_program_path = None):
        
        super().__init__(output_filename = output_filename,
                        Protein = Protein,
                        Ligand = Ligand,
                        protein_tpg_file = protein_tpg_file,
                        protein_prm_file = protein_prm_file,
                        MD_program_path = MD_program_path)
        #checking if I really have a protein and a ligand
        if Protein == None:
            raise ValueError("Portein can't be None type")

        elif Ligand == None:
            print("I found no organic ligands, I will optimize the protein on it's own")
        
        if self.output_filename == None:
            self.output_filename = f"{Protein.protein_id}_firstopt_orac.in"

        self.template = [
            "###############################################################",
            "#  Minimize Crystallographic structure form PDBank",
            "###############################################################",
            "",
            "#",
            "# Set MD cell and read pdb coordinates",
            "#",
            "&SETUP",
            
            self.write_box(), #"   CRYSTAL  150.0 150.0 150.0 90.0 90.0 90.0",

            f"   READ_PDB  {self.Protein.filename}",

            "&END",
            "",
            "#",
            "# read ASCII databases and build up solute",
            "#",
            "&PARAMETERS",

            f"   READ_TPG_ASCII {self.protein_tpg_file}",

            self.write_ligand_tpg_path(self.Ligand),

            f"   READ_PRM_ASCII {self.protein_prm_file}",

            self.write_ligand_prm_path(self.Ligand),

            "#TPGCYS",
            "ADD_TPG  SOLUTE  !! adds cys-cys",

            self.write_sulf_bond_string(self.Protein),

            "END",
            "   JOIN SOLUTE  !! primary structure",

            ' '+'\n '.join(Protein.seqres).lower(),

            self.get_ligand_name_from_tpg(self.Ligand),

            "   END",
            "&END",
            "",
            "#",
            "#  Simulation  Commands:  Minimize the structure",
            "#",
            "#",
            "&SIMULATION",
            "   MINIMIZE",
            "      CG  0.01",
            "   END",
            "&END",
            "",
            "#",
            "#  Cutoff for minimize is 7.0 A.",
            "#",
            "&POTENTIAL",
            "   UPDATE      20.0   1.5",
            "   CUTOFF  7.0",
            "   STRETCHING",
            "   QQ-FUDGE  0.83333",
            "   LJ-FUDGE  0.50",
            "&END",
            "",
            "#",
            "#  do 3 minimization step and intermediate printout every 5",
            "#",
            "&RUN",
            "   CONTROL      0",
            "   TIME         3.0",
            "   PRINT         5.0",
            "&END",
            "",
            "#overwriting the input pdb",
            f"# write final pdb file to {self.Protein.filename.rsplit('.', 1)[0].strip()}_optimized.pdb",
            "#",
            "&INOUT",
            f"   ASCII_OUTBOX    20.0 OPEN {self.Protein.filename.rsplit('.', 1)[0].strip()}_optimized.pdb",
            f"   PLOT FRAGMENT 1.0 OPEN {self.Protein.filename.rsplit('.', 1)[0].strip()}_optimized.xyz",
            "&END",
        ]


class OracSolvBoxInput(OracInput):
    """This is the class that contains the
    template to create the Orac input
    to insert the protein in a solvent box and make a MD
    simulation to reach aequilibrium at constant pressure"""

    #deactivating some 
    #BiopythonWarning
    import warnings
    import Bio.PDB.PDBExceptions
    warnings.simplefilter('ignore', Bio.PDB.PDBExceptions.PDBConstructionWarning)

    def __init__(self,
                input_filename = None,
                output_filename = None,
                Protein = None,
                Ligand = None,
                protein_tpg_file = None,
                protein_prm_file = None,
                MD_program_path = None,
                solvent_pdb = None):
        
        """takes the output_filename of the file on which to write
        and the istances of the present proteins and ligands"""

        super().__init__(input_filename = input_filename,
                        output_filename = output_filename,
                        Protein = Protein,
                        Ligand = Ligand,
                        protein_tpg_file = protein_tpg_file,
                        protein_prm_file = protein_prm_file,
                        MD_program_path= MD_program_path,
                        solvent_pdb = solvent_pdb)
        
        if output_filename == None:
            self.output_filename = f'{Protein.protein_id}_orac_solvbox.in'

        self.template = [
            "#&T NTHREADS    8   CACHELINE   16",
            "#&T NT-LEVEL1   2   CACHELINE   16",
            "#&T NT-LEVEL2   4   CACHELINE   16",
            "###############################################################",
            "#  Minimize Crystallographic structure from PDBank",
            "###############################################################",
            "",
            "! this is a comment",
            "!! two exclamation points: system-dependent section"
            "! one exclamation point: system indipendent section (same for all inputs)",
            "#",
            "# Set MD cell and read pdb coordinates",
            "#",
            "&SETUP",

            self.write_box(),

            "&END",
            "#",
            "# reads the force fields",
            "#",
            "&PARAMETERS",
            f"   READ_TPG_ASCII {self.protein_tpg_file} ! protein",

            self.write_ligand_tpg_path(self.Ligand),

            f"   READ_PRM_ASCII {self.protein_prm_file} ! protein",

            self.write_ligand_prm_path(self.Ligand),

            "#TPGCYS",
            "ADD_TPG  SOLUTE  !! adds cys-cys",

            self.write_sulf_bond_string(self.Protein),

            "END",
            "   JOIN SOLUTE  !! defines primary structure",

            ' '+'\n '.join(Protein.seqres).lower(),

            self.get_ligand_name_from_tpg(self.Ligand),

            "   END",

            f"   WRITE_TPGPRM_BIN  {Protein.protein_id}.tpgprm",

            "   JOIN SOLVENT   ! solvente",
            "       tip3",
            "   END",
            "&END",
            "&SOLUTE  !reads complex pdb file",

            f"   COORDINATES {self.Protein.filename}",

            "&END",
            "&SOLVENT    !! generates solvent grid",

            self.write_solvent_grid(),

            "    CELL  SC",
            "   INSERT 0.7",

            f"   COORDINATES {self.solvent_pdb}",

            "&END",
            "&SIMULATION  ! simulation parameters (same for all)",
            "   MDSIM",
            "   TEMPERATURE   280.0 20.0",
            "   ISOSTRESS PRESS-EXT 0.1 BARO-MASS 30.0",
            "   THERMOS",
            "      solute  10.0",
            "      solvent 10.0",
            "      cofm    10.0",
            "      temp_limit 1000.0",
            "   END",
            "&END",
            "&INTEGRATOR     ! integration parameters (same for all)",
            "   TIMESTEP       9.0",
            "   MTS_RESPA",
            "      step intra 2",
            "      step intra 2",
            "      step nonbond 2  5.1",
            "      step nonbond 5  7.8   reciprocal",
            "      step nonbond 1  10.0",
            "      test_times OPEN  G0.tt 20",
            "      very_cold_start 0.1",
            "  END",
            "&END",
            "&POTENTIAL  !! potential parameters",

            self.write_EWALD_PME(),

            self.write_ADD_STR_COM(),

            "   UPDATE      60.0   1.8",

            self.write_LINKED_CELL(),

            "   STRETCHING HEAVY",
            "   QQ-FUDGE  0.83333",
            "   LJ-FUDGE  0.50",
            "&END",
            "&RUN  ! run lenght (same for all)",
            "   CONTROL      0",
            "   PROPERTY     20000.0",
            "   REJECT       20000.0",
            "   TIME         10000.0",
            "   STEER        0.0 30000.0",
            "   PRINT        300.0",
            "&END",
            "",
            "#",
            "# write restart file every 60.0 (approximately)",
            "#",
            "&INOUT ! files I/O",
            "   RESTART",

            f"      write  15000.0  OPEN  {self.Protein.protein_id}_solventbox.rst",

            "   END",

            f"   ASCII   3000.0 OPEN {self.Protein.protein_id}_solventbox.pdb",

            f"   PLOT STEER_ANALYTIC  500.0  OPEN {self.Protein.protein_id}_solventbox.dat",

            "&END"                                              
        ]

    def write_solvent_grid(self):
        """Writes the informations about the solvent grid"""

        lx, ly, lz = self.orient.create_box(self.Protein.structure)

        nx, ny, nz = self.orient.create_solvent_grid(lx, ly, lz)

        string = f"   GENERATE   {nx}  {ny}  {nz}  !! the grid depends on the BOX"

        return string

    def write_EWALD_PME(self):
        """Writes the informations about the grid
        in the reciprocal reticle"""
        
        lx, ly, lz = self.orient.create_box(self.Protein.structure)

        pme_x, pme_y, pme_z = self.orient.create_recipr_solvent_grid(lx, ly, lz)

        string = f"   EWALD PME 0.37   {pme_x}  {pme_y}  {pme_z}   4 !! grid on the reciprocal"

        return string

    def write_ADD_STR_COM(self):
        """Writes the informations about the bounding between the protein and the ligands"""

        if self.Ligand == None:
            return ''

        #get the separeted structures of the protein and the ligands
        self.Protein.structure, ligand_structures = self.orient.separate_protein_ligand()

        #Make sure to save the structures of each ligand
        #may come in handy
        self.Ligand = self.orient.Ligand
        
        #Get the atom number range of each segment
        protein_atoms, ligand_atoms = self.orient.atom_numbers()

        string = []

        #get the distance of the centers of mass of the protein and the ligands
        #and create the string
        for i, ligand in enumerate(pipeline_functions.get_iterable(ligand_structures)):

            COM_Protein, COM_ligand, distance = self.orient.center_mass_distance(self.Protein.structure, ligand)

            tmp_string = ["   ADD_STR_COM   !! linker COM-COM legand protein",
                        f"       ligand     {ligand_atoms[i][1]}      {ligand_atoms[i][0]}",
                        f"       target    {protein_atoms[1]}      {protein_atoms[0]}",
                        f"       force   0.15     {distance}    {distance+0.0001}",
                        "   END"]
            
            tmp_string = '\n'.join(tmp_string)
            string.append(tmp_string)

        if len(string) > 1:
            string = '\n'.join(string)
        else:
            string = string[0]
        
        return string

    def write_LINKED_CELL(self):
        """writes the LINKED CELL string"""

        lx, ly, lz = self.orient.create_box(self.Protein.structure)

        string = f"   LINKED_CELL   {int(lx/3.0)}  {int(ly/3.0)}  {int(lz/3.0)}"

        return string










    



