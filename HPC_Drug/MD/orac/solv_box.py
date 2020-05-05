######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

from HPC_Drug.MD.orac import orac_input

import os

#deactivating some
#BiopythonWarning
import warnings
import Bio.PDB.PDBExceptions
warnings.simplefilter('ignore', Bio.PDB.PDBExceptions.PDBConstructionWarning)

class OracFirstOptimization(orac_input.OracInput):
    """
    Takes a HPC_Drug.structures.protein.Protein instance 
    with a self.pdb_file containing both the protein and the ligands (if any)
    the file MUST BE A PDB!

    Creates a solvent box and optimizes it NPT (constant pressure and temperature)

    The only pubblic methods are the contructor and execute()
    execute returns a HPC_Drug.structires.protein.Protein instance 
    """

    def __init__(self,
                Protein,
                solvent_pdb = None,
                MD_program_path = 'orac'):
        
        super().__init__(Protein = Protein,
                        solvent_pdb = solvent_pdb,
                        MD_program_path = MD_program_path)

        self.orac_in_file = os.getcwd() + f"/{self.Protein.protein_id}_solvbox_orac.in"
    
        self.output_pdb_file = os.getcwd() + f"/{self.Protein.protein_id}_solvbox.pdb"


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

            self._write_box(),

            "&END",
            "#",
            "# reads the force fields",
            "#",
            "&PARAMETERS",
            f"   READ_TPG_ASCII {self.Protein.tpg_file} ! protein",

            self._write_ligand_tpg_path(),

            f"   READ_PRM_ASCII {self.Protein.prm_file} ! protein",

            self._write_ligand_prm_path(),

            "#TPGCYS",
            "ADD_TPG  SOLUTE  !! adds cys-cys",

            self._write_sulf_bond_string(),

            "END",
            "   JOIN SOLUTE  !! defines primary structure",

            ' '+'\n '.join(Protein.seqres).lower(),

            self._get_ligand_name_from_tpg(),

            "   END",

            f"   WRITE_TPGPRM_BIN  {Protein.protein_id}.tpgprm",

            "   JOIN SOLVENT   ! solvent",
            "       tip3",
            "   END",
            "&END",
            "&SOLUTE  !reads complex pdb file",

            f"   COORDINATES {self.Protein.pdb_file}",

            "&END",
            "&SOLVENT    !! generates solvent grid",

            self._write_solvent_grid(),

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

            self._write_EWALD_PME(),

            self._write_ADD_STR_COM(),

            "   UPDATE      60.0   1.8",

            self._write_LINKED_CELL(),

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

            f"      write  15000.0  OPEN  {self.output_pdb_file.rsplit('.', 1)[0].strip()}.rst",

            "   END",

            f"   ASCII   3000.0 OPEN {self.output_pdb_file}",

            f"   PLOT STEER_ANALYTIC  500.0  OPEN {self.output_pdb_file.rsplit('.', 1)[0].strip()}.dat",

            "&END"                                              
        ]


