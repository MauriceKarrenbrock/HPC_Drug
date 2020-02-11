#containing functions and classes for the creation of a gromacs input
import Bio.PDB
import subprocess

from HPC_Drug import structures
from HPC_Drug import file_manipulation
from HPC_Drug import important_lists
from HPC_Drug import funcs4orac

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

                #give the right name to histidines in order to get the right protonation
                elif residue.resname in important_lists.hist_resnames:
                    
                    if ph < 6.0:
                        residue.resname = 'HIS'
                    
                    else:
                        residue.resname = 'HID'

                #if a cysteine is in a disulfide bond it is called CYX
                elif residue.resname in important_lists.cyst_resnames:
                    residue.resname = _disulf_cysteine_gromacs_substitutions(
                                                    input_list = Protein.substitutions_dict[res_id],
                                                    resname = residue.resname
                                                    )

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
                ligand_itp = None,
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
        
        self.ligand_itp = ligand_itp

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
    

    def execute(self,
                template = None,
                filename = None,
                output_file = None,
                MD_program_path = None,
                command_string = None):

        #must be adapted better
        if MD_program_path == None:
            MD_program_path = self.MD_program_path

        if template == None:
            template = self.template

        if filename == None:
            filename = self.output_filename
        
        if output_file == None:
            output_file = self.Protein.filename.rsplit('.', 1)[0].strip() + '_unnamed_output'
        
        if command_string == None:
            command_string = self.command_string

        filename = self.write_template_on_file(template, filename)

        print("Running Gromacs")
        r = subprocess.run([f'{MD_program_path} {command_string}'],
                        shell = True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        universal_newlines=True)

        print(r.stdout)
        print(r.stderr)

        if r.returncode != 0:
            raise Exception(f"Gromacs failure\n{r.stdout}\n{r.stderr}")
 
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
                input_filename = input_filename,
                output_filename = output_filename,
                Protein = Protein,
                Ligand = Ligand,
                protein_tpg_file = protein_tpg_file,
                solvent_model = solvent_model,
                MD_program_path = MD_program_path)
        #protein_tpg_file is actually the protein model chosen through choices_file
        
        self.command_string = f"{MD_program_path} pdb2gmx -f {Protein.filename} -o {Protein.protein_id}.gro -p {Protein.protein_id}.top < {output_filename}"

    
        self.template = [
                        f"{self.protein_tpg_file}",
                        f"{self.solvent_model}"]

    def execute(self,
            template = None,
            filename = None,
            output_file = None,
            MD_program_path = None,
            command_string = None,
            Protein = None):

        """Takes a protein instance and returns one"""

        choices_file = super().execute(template = template,
                        filename = filename,
                        output_file = output_file,
                        MD_program_path = MD_program_path,
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

        
        #Create the protein's gro and top file from it's pdb
        protein_gro_top = GromacsMakeProteinGroTop(input_filename = self.input_filename,
                                            output_filename = self.output_filename,
                                            Protein = self.Protein,
                                            Ligand = self.Ligand,
                                            protein_tpg_file = self.protein_tpg_file,
                                            solvent_model = self.solvent_model,
                                            MD_program_path = self.MD_program_path)

        Protein = protein_gro_top.execute()

        #join the protein and ligand pdb
        Protein = funcs4orac.join_ligand_and_protein_pdb(Protein = self.Protein,
                                                        Ligand = self.Ligand,
                                                        output_filename = self.output_filename)
        
        #make a new joined gro file
        r = subprocess.run([f'{MD_program_path} editconf -f {Protein.filename} -o {Protein.protein_id}_joined.gro'],
                        shell = True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        universal_newlines=True)

        print(r.stdout)
        print(r.stderr)

        if r.returncode != 0:
            raise Exception(f"Gromacs failure\n{r.stdout}\n{r.stderr}")

        Protein.gro_file = f"{Protein.protein_id}_joined.gro"

        #edit and rename the top file
        Protein = self._edit_top_file(Protein = Protein, Ligand = Ligand)
        
        return Protein

    def _edit_top_file(self, Protein, Ligand):
        """Adds the needed #include and other informations
        to the protein top file in order
        to include the ligand"""

        with open(Protein.top_file, 'r') as f:
            top = f.readlines()

        itp_insertion_string = f'#include "{Ligand.itp_file}"'
        compound_string = f'{Ligand.resname}                 1'
        with open(f'{Protein.protein_id}_joined.top', 'w') as f:
            top.insert(0, itp_insertion_string)
            top.insert(-1, compound_string)

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
                ligand_itp = None,
                solvent_model = None,
                MD_program_path = 'gmx'):

        super().__init__(input_filename = input_filename,
                        output_filename = output_filename,
                        Protein = Protein,
                        Ligand = Ligand,
                        protein_tpg_file = protein_tpg_file,
                        ligand_itp = ligand_itp,
                        solvent_model = solvent_model,
                        MD_program_path = MD_program_path)

        if self.output_filename == None:
            self.output_filename = f'{self.Protein.protein_id}_joined_optimized.gro'

        if self.input_filename == None:
            self.input_filename = f"{self.Protein.protein_id}_first_opt.mdp"

        #creates the .tpr and then optimizes the structure
        self.command_string = f"{self.MD_program_path} grompp -f {self.input_filename}.mdp -c {Protein.gro_file} -p {Protein.top_file} -o {Protein.protein_id}_joined_optimized.tpr -maxwarn 100 && gmx mdrun -s {Protein.protein_id}_joined_optimized.tpr  -c {self.output_filename}"

        self.template = ["integrator	= steep",
                        "nsteps		= 1000",
                        "emtol		= 100",
                        "emstep		= 0.01",
                        "nstxout 	= 1",
                        "nstenergy	= 1",
                        "rlist		= 1.0",
                        "coulombtype	= pme",
                        "vdw-type	= cut-off",
                        "rvdw		= 1.0",
                        "rcoulomb	= 1.0",
                        "constraints	= none",
                        "integrator	= steep",
                        "nsteps		= 1000",
                        "emtol		= 100",
                        "emstep		= 0.01",
                        "nstxout 	= 1",
                        "nstenergy	= 1",
                        "rlist		= 1.0",
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
                ligand_itp = None,
                solvent_model = None,
                MD_program_path = 'gmx'):

        super().__init__(input_filename = input_filename,
                        output_filename = output_filename,
                        Protein = Protein,
                        Ligand = Ligand,
                        protein_tpg_file = protein_tpg_file,
                        ligand_itp = ligand_itp,
                        solvent_model = solvent_model,
                        MD_program_path = MD_program_path)

        if self.output_filename == None:
            self.output_filename = f'{self.Protein.protein_id}_solv_box.gro'

        if self.input_filename == None:
            self.input_filename = f"{self.Protein.protein_id}_solv_box.mdp"

        #creates the .tpr and then optimizes the structure
        self.command_string = f"{self.MD_program_path} grompp -f {self.input_filename}.mdp -c {Protein.gro_file} -p {Protein.top_file} -o {Protein.protein_id}_joined_optimized.tpr -maxwarn 100 && gmx mdrun -s {Protein.protein_id}_joined_optimized.tpr  -c {self.output_filename}"

        self.template = []
   