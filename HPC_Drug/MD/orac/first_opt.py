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

class OracFirstOptimization(orac_input.OracInput):
    """
    Takes a HPC_Drug.structures.protein.Protein instance 
    with a self.pdb_file containing both the protein and the ligands (if any)
    the file MUST BE A PDB!

    Makes a fast first optimization

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

        self.orac_in_file = os.getcwd() + f"/{self.Protein.protein_id}_firstopt_orac.in"
    
        self.output_pdb_file = os.getcwd() + f"/{self.Protein.protein_id}_optimized.pdb"


        self.template = [
            "###############################################################",
            "#  Minimize Crystallographic structure form PDBank",
            "###############################################################",
            "",
            "#",
            "# Set MD cell and read pdb coordinates",
            "#",
            "&SETUP",
            
            self._write_box(), #"   CRYSTAL  150.0 150.0 150.0 90.0 90.0 90.0",

            f"   READ_PDB  {self.Protein.pdb_file}",

            "&END",
            "",
            "#",
            "# read ASCII databases and build up solute",
            "#",
            "&PARAMETERS",

            f"   READ_TPG_ASCII {self.Protein.tpg_file}",

            self._write_ligand_tpg_path(),

            f"   READ_PRM_ASCII {self.Protein.prm_file}",

            self._write_ligand_prm_path(),

            "#TPGCYS",
            "ADD_TPG  SOLUTE  !! adds cys-cys",

            self._write_sulf_bond_string(),

            "END",
            "   JOIN SOLUTE  !! primary structure",

            ' '+'\n '.join(Protein.seqres).lower(),

            self._get_ligand_name_from_tpg(),

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
            f"# write final pdb file to {self.output_pdb_file}",
            "#",
            "&INOUT",
            f"   ASCII_OUTBOX    20.0 OPEN {self.output_pdb_file}",
            f"   PLOT FRAGMENT 1.0 OPEN {self.output_pdb_file.rsplit('.', 1)[0].strip()}.xyz",
            "&END",
        ]

