######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
Hamiltonian Replica Exchange HREM for Orac
"""

import os
import math
import shutil

from HPC_Drug.MD.orac import orac_input
from HPC_Drug import important_lists
from HPC_Drug import funcs4slurm
from HPC_Drug import funcs4pbs

class HREMOracInput(orac_input.OracInput):
    """
    Prepares a working input to execute an Hamiltonian Replica Exchange on an
    HPC cluster with Orac
    """

    def __init__(self,
                input_filename = None,
                output_filename = None,
                Protein = None,
                Ligand = None,
                protein_tpg_file = None,
                protein_prm_file = None,
                MD_program_path = None,
                solvent_pdb = None,
                kind_of_processor = 'skylake',
                number_of_cores_per_node = 64):

        super().__init__(input_filename = input_filename,
                        output_filename = output_filename,
                        Protein = Protein,
                        Ligand = Ligand,
                        protein_tpg_file = protein_tpg_file,
                        protein_prm_file = protein_prm_file,
                        MD_program_path= MD_program_path,
                        solvent_pdb = solvent_pdb)
        
        self.orac_in_file = os.getcwd() + f'{Protein.protein_id}_HREM.in'

        self.kind_of_processor = kind_of_processor
        self.number_of_cores_per_node = number_of_cores_per_node

        self.template = [
            f"!!! THIS INPUT IS FOR {self.kind_of_processor} ARCHITECTURES",
            "#&T NTHREADS    8   CACHELINE   16",
            "#&T NT-LEVEL1   4   CACHELINE   16",
            "#&T NT-LEVEL2   4   CACHELINE   16",
            "&REM",

            self._write_BATTERIES_string(),
            
            "SETUP    1.0           0.2         1.0         1",
            "STEP 15.",
            "PRINT_DIAGNOSTIC 120.",
            "SEGMENT",

            self._write_SEGMENT_string(),
            
            "kind intra",
            "END",
            "PRINT 12000",
            
            f"PRINT_ENERGY 120.0 OPEN {self.Protein.protein_id}.rem",
            
            "&END",
            "&SETUP",

            self._write_box(),

            f"READ_PDB {self.Protein.pdb_file.rsplit('/', 1)[-1].strip()}",

            "&END",
            "&PARAMETERS",
            f"   READ_TPG_ASCII {self.Protein.tpg_file.rsplit('/', 1)[-1].strip()} ! protein",

            self._write_ligand_tpg_path(),

            f"   READ_PRM_ASCII {self.Protein.prm_file.rsplit('/', 1)[-1].strip()} ! protein",

            self._write_ligand_prm_path(),

            "#TPGCYS",
            "JOIN SOLUTE",

            ' '+'\n '.join(Protein.seqres).lower(),

            self._get_ligand_name_from_tpg(),
            
            "END",
            "JOIN SOLVENT",
            "tip3",
            "END",
            "&END",
            "&SOLVENT",
            "ADD_UNITS 6219",
            "&END",
            "&SIMULATION",
            "MDSIM",
            "TEMPERATURE   300.0 20.0",
            "ISOSTRESS PRESS-EXT 0.1 BARO-MASS 60.0",
            "THERMOS",
            "solute  30.0",
            "solvent 30.0",
            "cofm    30.0",
            "temp_limit 4000.0",
            "END",
            "&END",
            "&INTEGRATOR",
            "TIMESTEP       15.0",
            "MTS_RESPA",
            "step intra 2",
            "step intra 2",
            "step nonbond 2  5.1",
            "step nonbond 5  8.0   reciprocal",
            "step nonbond 1  10.0",
            "END",
            "&END",
            "&POTENTIAL",

            self._write_EWALD_PME(),

            self._write_ADD_STR_COM(),

            "END",
            "UPDATE      60.0   1.8",

            self._write_LINKED_CELL(),

            "STRETCHING HEAVY",
            "QQ-FUDGE  0.83333",
            "LJ-FUDGE  0.50",
            "&END",
            "&RUN",
            "CONTROL      0",
            "PROPERTY     20000000.0",
            "REJECT       30000.0",

            self._write_TIME_string(),

            "PRINT        1200.0",
            "&END",
            "",
            "#",
            "# write restart file every 60.0 (approximately)",
            "#",
            "&INOUT",
            "RESTART",
            
            f"write  45000.0 SAVE_ALL_FILES ../RESTART/{self.Protein.protein_id}",
            
            "END",
            
            f"PLOT FRAGMENT 9000.0 OPEN {self.Protein.protein_id}_rem.xyz",
            
            f"ASCII   150000.0 OPEN {self.Protein.protein_id}_rem.pdb",
            
            "&END",
            "&PROPERTIES",

            self._write_DEF_FRAGMENT_string(),

            "&END"
        ]

    def _write_ligand_tpg_path(self):
        """
        Private
        
        Writes the tpg file for any given ligand
        overwrites the superclass function
        """

        Ligand = self.Protein.get_ligand_list()
    
        if Ligand == None or Ligand == []:
            return ''

        tmp = []
        for ligand in Ligand:

            string = f"   READ_TPG_ASCII {ligand.tpg_file.rsplit('/', 1)[-1].strip()} !! ligand"

            tmp.append(string)

        string = '\n'.join(tmp)

        return string

    def _write_ligand_prm_path(self):
        """
        Private
        
        Writes the prm file for any given ligand
        overwrites the superclass function
        """

        Ligand = self.Protein.get_ligand_list()
    
        if Ligand == None or Ligand == []:
            return ''

        tmp = []
        for ligand in Ligand:

            string = f"   READ_PRM_ASCII {ligand.prm_file.rsplit('/', 1)[-1].strip()}  !! ligand"

            tmp.append(string)

        string = '\n'.join(tmp)

        return string

    def _write_SEGMENT_string(self):
        """
        Private
        """

        Ligand = self.Protein.get_ligand_list()
    
        if Ligand == None or Ligand == []:
            return ''

        Protein = self.Protein

        hot_residues = self.orient.get_hot_residues_for_rem(Protein = Protein, Ligand = Ligand, cutoff = 4.5, residue_dist = 10.0)

        DUMMY, ligand_atoms = self.orient.protein_ligand_atom_numbers(Protein = Protein, Ligand = Ligand)

        if len(ligand_atoms) > 0:
            for item in ligand_atoms:
                ligand_string = f"define {item[1]} {item[0]}    ! ligand"

        residue_string = ''
        for residue in hot_residues:
            if residue[0].strip().upper() != 'TIP':
                max_min = self.orient.atom_numbers(Protein = Protein, residue_id = residue[1])
                residue_string = residue_string + f"\ndefine {max_min[1]} {max_min[0]}    ! {residue[0]} {residue[1]}"

        SEGMENT_string = f"{ligand_string}\n{residue_string}"

        return SEGMENT_string

    def _write_DEF_FRAGMENT_string(self):
        """
        private

        writes the first and last atom of the protein ligand complex
        """

        Ligand = self.Protein.get_ligand_list()

        Protein = self.Protein

        max_min_atom = self.orient.get_first_last_atom_structure(Protein = Protein, Ligand = Ligand)

        return f"DEF_FRAGMENT   {max_min_atom[1]} {max_min_atom[0]}"

    def _get_ns_per_day(self):
        """
        Private
        
        Get's the number of ns per day that a kind_of_processor
        processor can process on the sistem
        """

        Ligand = self.Protein.get_ligand_list()

        Protein = self.Protein

        kind_of_processor = self.kind_of_processor

        number_of_atoms = self.orient.get_first_last_atom_structure(Protein = Protein, Ligand = Ligand)
        number_of_atoms = number_of_atoms[0]

        ns_per_day = ( 15000. / number_of_atoms ) * important_lists.processor_kind_ns_per_day_15000_atoms_for_cpu_only_runs[kind_of_processor]

        return ns_per_day

    def _write_TIME_string(self):
        """
        Private
        
        Get's the number of ns per day that a kind_of_processor
        processor can process on the sistem
        """

        time = self._get_ns_per_day()
        time = time * 1.E+6

        return f"TIME         {time}"

    def _get_BATTERIES(self):
        """
        Private
        
        Get's the number of batteries for REM
        """

        ns_per_day = self._get_ns_per_day()

        BATTERIES = math.ceil( 32. / ns_per_day )

        return BATTERIES

    def _write_BATTERIES_string(self):
        """
        Private
        """

        BATTERIES = self._get_BATTERIES()

        return f"BATTERIES    {BATTERIES}"

    def _write_workloadmanager_inputs(self):
        """private method called by self.execute()"""
        
        file_list = []
 
        slurm = funcs4slurm.SlurmInput(MD_input_file = self.orac_in_file,
                                    slurm_input_file = f'{self.orac_in_file.rsplit(".", 1)[0]}.slr',
                                    MD_program = 'orac',
                                    MD_calculation_type = 'rem',
                                    number_of_cores_per_node = self.number_of_cores_per_node,
                                    max_time = "24:00:00",
                                    ntasks = self._get_BATTERIES() * 8,
                                    cpus_per_task = 8,
                                    std_out = f'{self.orac_in_file.rsplit(".", 1)[0]}.out',
                                    std_err = f'{self.orac_in_file.rsplit(".", 1)[0]}.err')

        file_list.append(slurm.write())

        pbs = funcs4pbs.SlurmInput(MD_input_file = self.orac_in_file,
                                    slurm_input_file = f'{self.orac_in_file.rsplit(".", 1)[0]}.pbs',
                                    MD_program = 'orac',
                                    MD_calculation_type = 'rem',
                                    number_of_cores_per_node = self.number_of_cores_per_node,
                                    max_time = "24:00:00",
                                    ntasks = self._get_BATTERIES() * 8,
                                    cpus_per_task = 8,
                                    std_out = f'{self.orac_in_file.rsplit(".", 1)[0]}.out',
                                    std_err = f'{self.orac_in_file.rsplit(".", 1)[0]}.err')

        file_list.append(pbs.write())

        return file_list



    def execute(self, template = None, filename = None):
        """
        This execute method is different from the other ones
        because it doesn't run Orac, but writes an orac input to be used on a HPC cluster
        creates the needed directories and copies the needed files.
        and writes the workload manager input too (like SLURM)
        
        returns the absolute path of the orac input file
        """

        hrem_dir = f"{self.Protein.protein_id}_HREM"

        #creates the REM directory that will be copied to the HPC cluster
        os.makedirs(f"{hrem_dir}/RESTART", exist_ok=True)

        #for any existing ligand
        for ligand in self.Protein.get_ligand_list():
            #copy the ligand topology file to the new directory
            shutil.copy(os.path.abspath(os.path.expanduser(os.path.expandvars(ligand.tpg_file))), hrem_dir)
            #copy the ligand parameter file to the new directory
            shutil.copy(os.path.abspath(os.path.expanduser(os.path.expandvars(ligand.prm_file))), hrem_dir)

        #copy the protein topology file to the new directory
        shutil.copy(os.path.abspath(os.path.expanduser(os.path.expandvars(self.Protein.tpg_file))), hrem_dir)
        #copy the protein parameter file to the new directory
        shutil.copy(os.path.abspath(os.path.expanduser(os.path.expandvars(self.Protein.prm_file))), hrem_dir)

        #copy the protein pdb file to the new directory
        shutil.copy(os.path.abspath(os.path.expanduser(os.path.expandvars(self.Protein.pdb_file))), hrem_dir)


        #writes the orac input
        self._write_template_on_file()
        
        #copy the orac input file to the REM directory
        shutil.copy(os.path.abspath(os.path.expanduser(os.path.expandvars(self.orac_in_file))), hrem_dir)

        

        #write workload manager input for different workload managers (slurm pbs ...)
        workload_files = self._write_workloadmanager_inputs()

        #copy the workload manager files to the REM directory
        for item in workload_files:
            shutil.copy(os.path.abspath(os.path.expanduser(os.path.expandvars(item))), hrem_dir)


        return self.Protein


        

