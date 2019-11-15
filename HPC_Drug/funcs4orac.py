#contains all the classes and functions needed fot the orac MD program
import structures
import Bio.PDB
import file_manipulation
import important_lists
import pipeline_functions
import lib
import importlib_resources
import subprocess

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
                res_id = str(residue._id[1])

                #check if they are bounding a metal and are not a disulfide bond
                if res_id in Protein.substitutions_dict.keys():
                    if Protein.substitutions_dict[res_id][2] in important_lists.metals:

                        #make substitutions with the selected function
                        residue.resname = substitute(Protein.substitutions_dict[res_id])

                #give the right name to histidines in order to get the right protonation
                elif residue.resname in important_lists.hist_resnames:
                    
                    if ph < 5.0:
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
        resname = 'CYM'

    else:
        print(input_list[0], 'substitution not implemented, going on like nothing happened')
        resname = input_list[0]

    return resname

def custom_zinc_orac_substitutions(input_list):
    """It is a custom version of standard_gromacs_substitutions(list) function
    with custom names for zinc binding residues"""
    
    resname = standard_orac_substitutions(input_list)

    if resname == 'HID' and input_list[2] == 'ZN':
        resname = 'HDZ'
    
    elif resname == 'HIE' and input_list[2] == 'ZN':
        resname = 'HEZ'

    elif resname == 'CIM' and input_list[2] == 'ZN':
        resname = 'CIZ'

    return resname


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
        self.input_filename = input_filename
        if self.input_filename == None:
            #may call the mixing function in the future
            raise Exception('need an imput pdb')
        
        self.output_filename = output_filename

        self.protein_tpg_file = protein_tpg_file
        #if no path is given searches the standard amber tpg inside lib module
        if self.protein_tpg_file == None:
            with importlib_resources.path('lib', 'amber99sb-ildn.tpg') as path:
                self.protein_tpg_file = path
        
        self.protein_prm_file = protein_prm_file
        #if no path is given searches the standard amber prm inside lib module
        if self.protein_prm_file == None:
            with importlib_resources.path('lib', 'amber99sb-ildn.prm') as path:
                self.protein_tpg_file = path

        self.solvent_pdb = solvent_pdb
        #if no path is given searches the standard water.pdb inside lib module
        if self.solvent_pdb == None:
            with importlib_resources.path('lib', 'water.pdb') as path:
                self.protein_tpg_file = path
        
        self.MD_program_path = MD_program_path
        if self.MD_program_path == None:
            raise Exception('Need a MD_program_path (for example ~/ORAC/trunk/bin/orac)')

        self.template = None

    def write_sulf_bond_string(self, Protein = None):
        """Writes the part of the template inherent to the disulf bonds"""

        if Protein == None:
            Protein = self.Protein
        
        tmp_sulf = []
        
        for bond in Protein.sulf_bonds:

            string = f"   bond 1sg 2sg residue   {bond[0]}  {bond[1]}"
            tmp_sulf.append(string)
        
        string = '\n'.join(tmp_sulf)

        return string

    def write_ligand_tpg_path(self, Ligand = None):
        """Writes the tpg file for any given ligand"""

        if Ligand == None:
            Ligand = self.Ligand

        tmp = []
        for ligand in pipeline_functions.get_iterable(Ligand):

            string = f"   READ_TPG_ASCII {ligand.topology_file}                                !! ligando"

            tmp.append(string)

        string = '\n'.join(tmp)

        return string

    def write_ligand_prm_path(self, Ligand = None):
        """Writes the prm file for any given ligand"""

        import pipeline_functions

        if Ligand == None:
            Ligand = self.Ligand

        tmp = []
        for ligand in pipeline_functions.get_iterable(Ligand):

            string = f"   READ_PRM_ASCII {ligand.param_file}                                !! ligando"

            tmp.append(string)

        string = '\n'.join(tmp)

        return string

    def get_ligand_name_from_tpg(self, Ligand = None):
        """Gets the name of the ligand from the tpg file
        for any given ligand"""

        if Ligand == None:
            Ligand = self.Ligand

        residue_strings = []

        for ligand in pipeline_functions.get_iterable(Ligand):
            with open(ligand.topology_file, 'r') as f:

                for line in f:

                    if 'RESIDUE' in line:
                        
                        string = f"      {line.split()[1].strip()}    !! nome del ligando nel file tpg"

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
                f.write(line)

        return filename

    def execute(self, template = None, filename = None, pdb_file = None, MD_program_path = None):

        if MD_program_path == None:
            MD_program_path = self.MD_program_path

        filename = self.write_template_on_file(template, filename)

        subprocess.run(f'{MD_program_path} {filename}', shell = True, check = True)

        return pdb_file  





class OracFirstOptimization(OracInput):
    """Takes the input file for the first optimization
    of the protein + ligand function (the 2 or more pdb files must already be mixed)
    
    Takes a Protein and any ligand instance
    the ligand needs an already calculated .tpg and .prm file

    makes the first optimization"""

    def __init__(self,
                input_filename = None,
                output_filename = None,
                Protein = None,
                Ligand = None,
                protein_tpg_file = None,
                protein_prm_file = None):
        
        super().__init__(input_filename,
                        output_filename,
                        Protein,
                        Ligand,
                        protein_tpg_file,
                        protein_prm_file)
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
            "   CRYSTAL  150.0 150.0 150.0 90.0 90.0 90.0",

            f"   READ_PDB  {self.input_filename}",

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
            "ADD_TPG  SOLUTE  !! aggiunge cys-cys",

            self.write_sulf_bond_string(self.Protein),

            "END",
            "   JOIN SOLUTE  !! definisce struttura primaria",
            
            self.get_ligand_name_from_tpg(self.Ligand),

            ' '+'\n '.join(Protein.seqres).lower(),
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
            "#  do 20 minimization step and intermediate printout every 5",
            "#",
            "&RUN",
            "   CONTROL      0",
            "   TIME         10.0",
            "   PRINT         5.0",
            "&END",
            "",
            "#overwriting the input pdb",
            f"# write final pdb file to {self.input_filename}",
            "#",
            "&INOUT",
            f"   ASCII_OUTBOX    20.0 OPEN {self.input_filename}",
            f"   PLOT FRAGMENT 1.0 OPEN {self.input_filename}.xyz",
            "&END",
        ]


class OracSolvBoxInput(OracInput):
    """This is the class that contains the
    template to create the Orac input
    to insert the protein in a solvent box and make a MD
    simulation to reach aequilibrium at constant pressure"""

    #####################################
    #DA FINIRE
    #######################################

    def __init__(self,
                input_filename = None,
                output_filename = None,
                Protein = None,
                Ligand = None,
                protein_tpg_file = None,
                protein_prm_file = None,
                solvent_pdb = None):
        
        """takes the output_filename of the file on which to write
        and the istances of the present proteins and ligands"""
        
        #########################
        #DA FINIRE!!!!!!!!!!
        #####################

        super().__init__(input_filename,
                        output_filename,
                        Protein,
                        Ligand,
                        protein_tpg_file,
                        protein_prm_file,
                        solvent_pdb)
        
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
            "!N.B questo e' un commento",
            "!! due punti esclamativi indicano che  la sezione di",
            "!! riferimento e' system-dependent",
            "! Un punto esclamativo indica che la sezione e' uguale per tutti i",
            "! sistemi",
            "#",
            "# Set MD cell and read pdb coordinates",
            "#",
            "&SETUP",
            "   CRYSTAL    82.13    75.25    58.69  !! BOX di simulazione",
            "&END",
            "#",
            "# legge i force field",
            "#",
            "&PARAMETERS",
            f"   READ_TPG_ASCII {self.protein_tpg_file} ! proteina",

            self.write_ligand_tpg_path(self.Ligand),

            f"   READ_PRM_ASCII {self.protein_prm_file} ! proteina",

            self.write_ligand_prm_path(self.Ligand),

            "#TPGCYS",
            "ADD_TPG  SOLUTE  !! aggiunge cys-cys",

            self.write_sulf_bond_string(self.Protein),

            "END",
            "   JOIN SOLUTE  !! definisce struttura primaria",
            
            self.get_ligand_name_from_tpg(self.Ligand),

            ' '+'\n '.join(Protein.seqres).lower(),

            "   END",

            f"   WRITE_TPGPRM_BIN  {Protein.protein_id}.tpgprm",

            "   JOIN SOLVENT   ! solvente",
            "       tip3",
            "   END",
            "&END",
            "&SOLUTE  !legge le coordinate PDB del complesso.",
            f"   COORDINATES {self.input_filename}",
            "&END",
            "&SOLVENT    !! genera solvente su una griglia",
            "   GENERATE   25  23  18  !! la griglia dipende dal BOX",
            "    CELL  SC",
            "   INSERT 0.7",
            "   COORDINATES ../../lib/water.pdb",
            "&END",
            "&SIMULATION                  ! parametri di simulazione (uguali per tutti)",
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
            "&INTEGRATOR     ! parametri di integrazione  (uguali per tutti)",
            "   TIMESTEP       9.0",
            "   MTS_RESPA",
            "      step intra 2",
            "      step intra 2",
            "      step nonbond 2  5.1",
            "      step nonbond 5  7.8   reciprocal",
            "      step nonbond 1  10.0",
            "      test_times OPEN  G0.tt 20",
            "      energy_then_die",
            "  END",
            "&END",
            "&POTENTIAL  !! parametri del potenziale",
            "   EWALD PME 0.37   64  64  48   4 !! griglia sul reciproco",
            "   ADD_STR_COM   !! linker COM-COM legando proteina",
            "       ligand     1      84",
            "       target    85    5830",
            "       force   0.15     8.9221    8.9222",
            "   END",
            "   UPDATE      60.0   1.8",
            "   LINKED_CELL   27  25  20",
            "   STRETCHING HEAVY",
            "   QQ-FUDGE  0.83333",
            "   LJ-FUDGE  0.50",
            "&END",
            "&RUN  ! lunghezza  del run (uguali per tutti)",
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
            "&INOUT ! files I/O  (uguali per tutti)",
            "   RESTART",
            "      write  15000.0  OPEN  md.rst",
            "   END",
            "   ASCII   3000.0 OPEN md_1.pdb",
            "   PLOT STEER_ANALYTIC  500.0  OPEN com.dat",
            "&END"                                              
        ]

    



