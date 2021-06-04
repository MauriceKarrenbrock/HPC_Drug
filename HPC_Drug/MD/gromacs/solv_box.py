######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

import os
import shutil
try: # python>=3.7

    import importlib.resources as importlib_resources

except ImportError: # python<=3.6

    import importlib_resources

from HPC_Drug.MD.gromacs import gromacs_input
from HPC_Drug.files_IO import write_on_files
from HPC_Drug.files_IO import read_file


class GromacsSolvBoxInput(gromacs_input.GromacsInput): 
    """Creates and optimizes a water box around the protein
    
    Parameters
    --------------
    Protein : HPC_Drug.structures.protein.Protein
    MD_program_path : str, default= looks for gmx in the path and current directory
        it must be the absolute path to the executable
    box_borders : float, default=1.2
        how far the box wall should be from the most outer atom of the structure in NM
    box_shape : str, default=cubic
        the shape of the box, the only accepted options are cubic and triclinic (rectangular, 90 angles)
    """

    def __init__(self,
                Protein = None,
                MD_program_path = 'gmx',
                box_borders = 1.2,
                box_shape = 'cubic'):

        super().__init__(Protein = Protein,
                        MD_program_path = MD_program_path)

        self.output_gro_file = os.getcwd() + "/" + f"{self.Protein.protein_id}_solvbox.gro"
        
        self.output_pdb_file = os.getcwd() + "/" +  f"{self.Protein.protein_id}_solvbox.pdb"

        self.mdp_file = os.getcwd() + "/" + f"{self.Protein.protein_id}_solvbox.mdp"

        self.output_tpr_file = os.getcwd() +  "/" + f"{self.Protein.protein_id}_solvbox.tpr"
        
        if box_borders is None:

            box_borders = 1.2
        
        self.box_borders = box_borders

        if box_shape is None:

            box_shape = 'cubic'

        elif box_shape not in ('cubic', 'triclinic'):

            raise ValueError(f'box_shape can only be cubic or triclinic, not {box_shape}')

        self.box_shape = box_shape

        #creates the .tpr of the structure + the box of water and then optimizes the system
        self.command_string = [
            [f"{self.MD_program_path}", "editconf", "-f", f"{self.Protein.gro_file}", "-d", f"{self.box_borders}", "-bt", f"{box_shape}", "-angles", "90", "90", "90", "-o", f"{self.output_gro_file}"],
            [f"{self.MD_program_path}", "solvate", "-cp", f"{self.output_gro_file}", "-p", f"{self.Protein.top_file}", "-o", f"{self.output_gro_file}"],
            [f"{self.MD_program_path}", "grompp", "-f", f"{self.mdp_file}", "-c", f"{self.output_gro_file}", "-p", f"{self.Protein.top_file}", "-maxwarn", "100", "-o", f"{self.output_tpr_file}"],
            [f"{self.MD_program_path}", "mdrun", "-s", f"{self.output_tpr_file}", "-c", f"{self.output_gro_file}", '-ntmpi', '1']
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
                        "dt                       = 0.0005",
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
                        "nstxtcout                = 1000",
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
                        "pcoupl                   = Parrinello-Rahman",
                        "pcoupltype               = Isotropic",
                        "nstpcouple               = -1",
                        "; Time constant (ps), compressibility (1/bar) and reference P (bar)",
                        "tau-p                    = 1.0",
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
                        "constraints              = h-bonds",
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

    # def _edit_top_file(self):
    #     """
    #     private
        
    #     Adds the needed #include for water
    #     to the protein top file
    #     """

    #     top = read_file.read_file(file_name = self.Protein.top_file)

    #     itp_insertion_string = f'#include "{self.solvent_model}"'

    #     for i in range(len(top)):
    #         if top[i].strip() == "" and top[i].strip()[0] == ";":
    #             pass
    #         elif top[i].strip().replace(" ", "").split(";")[0] == '[moleculetype]':
    #             top[i] = itp_insertion_string + '\n' + top[i]
    #             break


    #     write_on_files.write_file(lines = top, file_name = self.Protein.top_file)

    # def _pre_run_hook(self):
    #     """
    #     private
        
    #     overwriting the pre run hook method
    #     """

    #     self._edit_top_file()




class OptimizeOnlyWaterBox(gromacs_input.GromacsInput):
    """
    Optimizes a box of solvent
    as deafault uses a copy of the one found in HPC_Drug.lib
    """

    def __init__(self,
                solvent_box = None,
                force_fileld = "amber99sb-ildn",
                solvent_model = "spce",
                MD_program_path = "gmx"):

        self.solvent_box = solvent_box
        if self.solvent_box is None:
            with importlib_resources.path('HPC_Drug.lib', 'only_water.gro') as path:
                shutil.copy(str(path.resolve()), os.getcwd())

            self.solvent_box = os.getcwd() + "/" + "only_water.gro"

        self.force_fileld = force_fileld

        self.solvent_model = solvent_model

        self.MD_program_path = MD_program_path

        self.output_top_file = os.getcwd() + "/" + f"only_water_{self.solvent_model}.top"

        self.output_gro_file = os.getcwd() + "/" + f"only_water_{self.solvent_model}.gro"

        self.mdp_file = os.getcwd() + "/" + f"only_water_{self.solvent_model}_opt.mdp"

        self.output_tpr_file = os.getcwd() + "/" + f"only_water_{self.solvent_model}_opt.tpr"

        self.command_string = [
            [f"{self.MD_program_path}", "grompp", "-f", f"{self.mdp_file}", "-c", f"{self.output_gro_file}", "-p", f"{self.output_top_file}", "-maxwarn", "100", "-o", f"{self.output_tpr_file}"],
            [f"{self.MD_program_path}", "mdrun", "-s", f"{self.output_tpr_file}", "-c", f"{self.output_gro_file}", '-ntmpi', '1']
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
                        "dt                       = 0.0005",
                        "nsteps                   = 60000",
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
                        "nstxtcout                = 1000",
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
                        "pcoupl                   = Parrinello-Rahman",
                        "pcoupltype               = Isotropic",
                        "nstpcouple               = -1",
                        "; Time constant (ps), compressibility (1/bar) and reference P (bar)",
                        "tau-p                    = 1.0",
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
                        "constraints              = h-bonds",
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


    def _create_solvent_top(self):

        command_string = [f"{self.MD_program_path}", "pdb2gmx",
            "-ff", f"{self.force_fileld}",
            "-water", f"{self.solvent_model}",
            "-f", f"{self.solvent_box}",
            "-o", f"{self.output_gro_file}",
            "-p", f"{self.output_top_file}"]

        self._interact_with_gromacs(string = command_string)


    def execute(self):
        """
        Returns the updated gro and top files

        return gro_file, top_file
        """

        self._write_template_on_file()

        self._create_solvent_top()

        self._interact_with_gromacs()

        return self.output_gro_file, self.output_top_file



        