######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

import os
import os.path
import subprocess
import copy
import shutil
import math

from HPC_Drug.auxiliary_functions import path
from HPC_Drug.files_IO import read_file, write_on_files
from HPC_Drug.MD import workload_managers
from HPC_Drug.MD.gromacs import gro2pdb
from HPC_Drug.PDB import biopython, prody
from HPC_Drug import orient


class FSDAMInputProteinLigand(object):

    def __init__(self,
                HREM_dir = os.getcwd(),
                delta_distance = 4.,
                gromacs_path = "gmx",
                cpus_per_node = 64,
                use_GPU = False,
                GPU_per_node = 1):

        #the root directory of the HREM
        #default is working directory
        self.HREM_dir = HREM_dir
        
        self.delta_distance = delta_distance

        self.gromacs_path = path.absolute_programpath(program = gromacs_path)

        self.cpus_per_node = cpus_per_node

        self.use_GPU = use_GPU

        self.GPU_per_node = GPU_per_node


    def _create_restart_configs(self, reference_distance, have_reference, fsdam_dir):
        """
        Private

        Makes the starting files for fs-dam and writes them in the RESTART dir

        But first checks if the ligand remained in place 
        """

        i = 0
        while os.path.exists(f"BATTERY{i}"):

            subprocess.run(f"{self.gromacs_path} distance -f BATTERY{i}/scaled0/traj_comp.xtc -s BATTERY{i}/scaled0/topol.tpr -select 'com of group \"Protein\" plus com of group \"S7V\"' -oh histogramm_BATTERY{i}.xvg -len 1 -tol 1.5 -binw 0.005 -xvg none -rmpbc yes -oav average_dist_BATTERY{i}.xvg",
                        shell = True)

            with open(f"average_dist_BATTERY{i}.xvg", "r") as f:

                lines = f.readlines()

            if not have_reference:

                reference_distance = float(lines[0].strip().split()[1].strip())

            #my measures are in Angstrom, gromacs measures in nm
            actual_distance = float(lines[-1].strip().split()[1].strip()) * 10.
            if abs( actual_distance - reference_distance ) > self.delta_distance:

                #if the ligand moved to much raises a RuntimeError
                raise RuntimeError(f"The ligand moved away from the docked position for more than {self.delta_distance} Angstrom")

            else:

                subprocess.run(f"echo non-Water | {self.gromacs_path} trjconv -f BATTERY{i}/scaled0/traj_comp.xtc -s BATTERY{i}/scaled0/topol.tpr -pbc mol -b 80 -o {fsdam_dir}/BATTERY{i}_.gro -sep yes -ur compact", shell = True)
                
                i = i + 1


    def _create_itp_files(self, input_itp_filename, ligand_annihilation_itp_filename):


        #make one itp file ligand annihilation
        untouched_itp = read_file.read_file(file_name = input_itp_filename)

        #add dummy atoms and create [ nonbond_params ]

        nonbond_params = [
            "[ nonbond_params ]\n",
            "; i  j    func  sigma           epsilon\n"
        ]

        atoms = []

        is_atomtype = False
        for i in range(len(untouched_itp)):

            if untouched_itp[i].strip().replace(" ", "") == "[atomtypes]" and is_atomtype == False:

                is_atomtype = True

            elif (is_atomtype == True) and (untouched_itp[i].strip()[0] == "["):

                #This will be the last thing to happen (after the next if)

                for i in range(len(atoms)):

                    for j in range(i, len(atoms)):

                        #DUM_atom i
                        tmp_nonbond_params = "DUM_" + atoms[i][0].strip()
                        #DUM_atom j
                        tmp_nonbond_params = tmp_nonbond_params + " " + "DUM_" + atoms[j][0].strip()
                        #func (always 1)
                        tmp_nonbond_params = tmp_nonbond_params + " " + "1"
                        #sigma (arithmetic awg of the two sigmas)
                        tmp_nonbond_params = tmp_nonbond_params + " " + f"{( float( atoms[i][5].strip() ) + float( atoms[j][5].strip() ) ) / 2.}"
                        #epsilon (geometric awg of the two epsilons)
                        tmp_nonbond_params = tmp_nonbond_params + " " + f"{( float( atoms[i][6].strip() ) * float( atoms[j][6].strip() ) ) ** 0.5}"
                        #newline
                        tmp_nonbond_params = tmp_nonbond_params + "\n"

                        #append to the cumulative ones
                        nonbond_params.append(tmp_nonbond_params)

                nonbond_params.append("\n\n")

                #add the nonbond_params in untouched_itp 
                untouched_itp[i:i] = nonbond_params

                break


            if is_atomtype and (untouched_itp[i].strip() != "") and (untouched_itp[i].strip()[0] != ";"):

                line = untouched_itp[i].strip().split()

                tmp_DUM = ["DUM_" + line[0].strip(), line[1].strip(), line[2].strip(), "0.0", line[4].strip(), "0.0", "0.0", "\n"]
                
                tmp_DUM = " ".join(tmp_DUM)

                #every atom will be daclared with its DUM_ version
                untouched_itp[i] = untouched_itp[i].strip() + "\n" + tmp_DUM

                #keeping the normal atoms in order to do the nonbond_params later
                atoms.append(line)


        #itp for ligand gets annihilated (protein ligand system)
        ligand_annihilation_itp = untouched_itp

        is_atoms = False
        for i in range(len(ligand_annihilation_itp)):

            if ligand_annihilation_itp[i].strip().replace(" ", "") == "[atoms]" and is_atoms == False:

                is_atoms = True

            elif is_atoms and ligand_annihilation_itp[i].strip()[0] == "[":

                break

            if is_atoms and (ligand_annihilation_itp[i].strip() != "") and (ligand_annihilation_itp[i].strip()[0] != ";"):
                
                tmp_line = ligand_annihilation_itp[i].strip().split()

                ligand_annihilation_itp[i] = ligand_annihilation_itp[i].strip()

                #add DUM_ state b
                ligand_annihilation_itp[i] = ligand_annihilation_itp[i] + " " + "DUM_" + tmp_line[1]

                ligand_annihilation_itp[i] = ligand_annihilation_itp[i] + " " + "0.0" + tmp_line[7]

                ligand_annihilation_itp[i] = ligand_annihilation_itp[i] + "\n"


        #writing it in the right place
        write_on_files.write_file(lines = ligand_annihilation_itp, file_name = ligand_annihilation_itp_filename)
        

    def _gro2pdb(self, gro, pdb):
        
        gro2pdb.gro2pdb(gro_file = gro, pdb_file = pdb, gromacs_path = self.gromacs_path)

    def _get_heavy_atom(self, ligand_pdb_file):
        """
        private

        given a pdb_file containing only the ligand finds the atom number of the nearest
        ligand atom to the center of mass (COM) that is not an hydrogen
        """

        orient_obj = orient.Orient()

        structure = biopython.parse_pdb("idid", ligand_pdb_file)
        
        ligand_com = orient_obj.center_of_mass(structure)

        #get the distance from the first atom
        atoms = structure.get_atoms()
        for atom in atoms:
            min_dist = ( atom.coord[0] - ligand_com.coord[0] ) ** 2
            min_dist = min_dist + ( atom.coord[1] - ligand_com.coord[1] ) ** 2
            min_dist = min_dist + ( atom.coord[2] - ligand_com.coord[2] ) ** 2
            min_dist = min_dist ** 0.5

            heavy_atom = atom.serial_number

            break

        #now check for the nearest atom that is not an Hydrogen
        H = ("H", "h")
        atoms = structure.get_atoms()
        for atom in atoms:

            if atom.name.strip()[0] not in H:

                distance = ( atom.coord[0] - ligand_com.coord[0] ) ** 2
                distance = distance + ( atom.coord[1] - ligand_com.coord[1] ) ** 2
                distance = distance + ( atom.coord[2] - ligand_com.coord[2] ) ** 2
                distance = distance ** 0.5

                if distance < min_dist:

                    heavy_atom = atom.serial_number

                    min_dist = distance

        return heavy_atom

    def _get_only_ligand_pdb(self, pdb_file, resnum = None, resname = None):
    
        struct = prody.parse_pdb(pdb_file)

        selecter = prody.ProdySelect(struct)

        if resnum is not None:

            struct = selecter.resnum(resnum)

        elif resname is not None:

            struct = selecter.resname(resname)

        else:
            raise ValueError("I need a resnum or a resname, cannot both be None")

        prody.write_pdb(structure = struct, file_name = "tmp_only_ligand.pdb")

        return "tmp_only_ligand.pdb"

    def _make_heavy_atom_heavy_in_itp(self, heavy_atom, itp_file):

        lines = read_file.read_file(file_name = itp_file)
        
        #auxiliary bool variable
        is_atoms = False
        for i in range(len(lines)):

            #if the line is empty i can go on with the for loop
            if lines[i].strip() != "":
                tmp_line = lines[i].strip().split()
            else:
                continue

            #check if we are in the atoms part
            if tmp_line[0].replace(" ", "") == "[atoms]":
                is_atoms = True

            #end of atoms part
            elif tmp_line[0][0] == "[":
                is_atoms = False

            if is_atoms:

                #if it is not a comment
                if tmp_line[0][0] != ";":

                    if int(tmp_line[0].strip()) == heavy_atom:

                        tmp_line[7] = tmp_line[7].strip() + "E+20"

                        lines[i] = " ".join(tmp_line)

                        break

        else:
            raise RuntimeError(f"Could not find the heavy atom {heavy_atom} in this itp file")

        write_on_files.write_file(lines = lines, file_name = itp_file) 


    def _set_velocity_to_zero_4_heavy_atom(self,
                                        gro_file,
                                        heavy_atom,
                                        ligand_resnum = None,
                                        ligand_resname = None):

        if ligand_resnum is None and ligand_resname is not None:        

            ligand_resname = (ligand_resname.upper(), ligand_resname.lower())

        elif ligand_resname is None:

            raise ValueError("I need a resnum or resname the cannot both be None")
        

        lines = read_file.read_file(file_name = gro_file)

        #check the first atom number of the ligand in the joined gro file

        for line in lines:

            if ligand_resnum is not None and int(line[:6].strip()) == ligand_resnum:

                first_ligand_atom = int(line[16:21].strip())
                
                break

            elif ligand_resname is not None and int(line[:6].strip()) in ligand_resname:

                first_ligand_atom = int(line[16:21].strip())
                
                break

        heavy_atom = heavy_atom + first_ligand_atom

        for i in range(len(lines)):

            if int(lines[i][16:21].strip()) == heavy_atom:

                lines[i] = lines[i][:44] + "{0:8.4f}{0:8.4f}{0:8.4f}\n".format(0.)

        write_on_files.write_file(lines = lines, file_name = gro_file)
        
    def _prepare_input(self, ligand_resname, fsdam_dir, top_file):
        
        prepare_input = PrepareInputProteinLigand(fsdam_dir = fsdam_dir,
                            ligand_resname = ligand_resname,
                            top_file = top_file,
                            cpus_per_node = self.cpus_per_node,
                            use_GPU = self.use_GPU,
                            GPU_per_node = self.GPU_per_node)

        prepare_input.execute()

    
    def execute(self):

        #move to the HREM directory
        if self.HREM_dir != os.getcwd():
            
            changed_dir = True

            old_dir = os.getcwd()

            os.chdir(self.HREM_dir)

        else:

            changed_dir = False

            old_dir = None


        fsdam_dir = "RESTART"

        os.makedirs(fsdam_dir, exist_ok=True)

        #lazily copy all the itp files in the new dir
        for i in os.listdir(os.getcwd()):

            if i[-4:] == ".itp":

                shutil.copy(i, fsdam_dir)

        #making some defaults
        useful_info = {
            "ligand_resname" : None,
            "ligand_resnum" : None,
            "top_file" : "topol.top",
            "ligand_itp" : "LIG.itp",
            "reference_distance" : None
        }

        lines = read_file.read_file(file_name = "important_info.dat")
        for line in lines:

            line = line.strip()
            line = line.split("=")
            useful_info[line[0].strip()] = line[1].strip()

        if useful_info["reference_distance"] is None:
            have_reference = False

        else:
            useful_info["reference_distance"] = float(useful_info["reference_distance"])
            have_reference = True

        if useful_info["ligand_resnum"] is not None:
            useful_info["ligand_resnum"] = int(useful_info["ligand_resnum"])

        self._create_restart_configs(reference_distance = useful_info["reference_distance"],
                                    have_reference = have_reference,
                                    fsdam_dir = fsdam_dir)

        
        #find the nearest not H atom to the center of mass
        #make it heavy ~ E+20 Da, and set its velocity to 0 in any gro file
        #inside fsdam_dir
        for i in os.listdir(fsdam_dir):

            if i[-4:] == ".gro":
                
                self._gro2pdb(gro = fsdam_dir + "/" + i, pdb = "tmp_protein_ligand.pdb")

                break


        only_ligand_pdb = self._get_only_ligand_pdb(pdb_file = "tmp_protein_ligand.pdb",
                                                resnum = useful_info["ligand_resnum"],
                                                resname = useful_info["ligand_resname"])

        
        heavy_atom = self._get_heavy_atom(only_ligand_pdb)

        self._make_heavy_atom_heavy_in_itp(heavy_atom = heavy_atom, itp_file = useful_info['ligand_itp'])

        for i in os.listdir(fsdam_dir):
        
            if i[-4:] == ".gro":
                self._set_velocity_to_zero_4_heavy_atom(gro_file = fsdam_dir + "/" + i,
                                                    heavy_atom = heavy_atom,
                                                    ligand_resnum = useful_info["ligand_resnum"],
                                                    ligand_resname = useful_info["ligand_resname"])


        #creates the two new itp files (ligand creation and annihilation)
        self._create_itp_files(input_itp_filename = useful_info['ligand_itp'],
                            ligand_annihilation_itp_filename = f"{fsdam_dir}/{useful_info['ligand_itp']}")


        #copy the protein ligand top in the protein_ligand dir
        shutil.copy(useful_info["top_file"], fsdam_dir)

        #prepares the input files for the 
        self._prepare_input(ligand_resname = useful_info["ligand_resname"],
                            fsdam_dir = fsdam_dir, top_file = useful_info["top_file"])


        #go back to the old working directory (attention, if something fails the working dir will remain in the
        # HREM_dir!! (use try: except: in case needed))

        if changed_dir:

            os.chdir(old_dir)


        



class FSDAMInputOnlyLigand(object):

    def __init__(self,
                HREM_dir = os.getcwd(),
                gromacs_path = "gmx",
                cpus_per_node = 64,
                use_GPU = False,
                GPU_per_node = 1):

        #the root directory of the HREM
        #default is working directory
        self.HREM_dir = HREM_dir
        
        self.gromacs_path = path.absolute_programpath(program = gromacs_path)

        self.cpus_per_node = cpus_per_node

        self.use_GPU = use_GPU

        self.GPU_per_node = GPU_per_node


    def _create_itp_files(self, input_itp_filename, itp_filename):


        #make the two itp files, ligand annihilation and ligand creation
        untouched_itp = read_file.read_file(file_name = input_itp_filename)

        #add dummy atoms and create [ nonbond_params ]
        nonbond_params = [
            "[ nonbond_params ]\n",
            "; i  j    func  sigma           epsilon\n"
        ]

        atoms = []

        is_atomtype = False
        for i in range(len(untouched_itp)):

            if untouched_itp[i].strip().replace(" ", "") == "[atomtypes]" and is_atomtype == False:

                is_atomtype = True

            elif (is_atomtype == True) and (untouched_itp[i].strip()[0] == "["):

                #This will be the last thing to happen (after the next if)

                for i in range(len(atoms)):

                    for j in range(i, len(atoms)):

                        #DUM_atom i
                        tmp_nonbond_params = "DUM_" + atoms[i][0].strip()
                        #DUM_atom j
                        tmp_nonbond_params = tmp_nonbond_params + " " + "DUM_" + atoms[j][0].strip()
                        #func (always 1)
                        tmp_nonbond_params = tmp_nonbond_params + " " + "1"
                        #sigma (arithmetic awg of the two sigmas)
                        tmp_nonbond_params = tmp_nonbond_params + " " + f"{( float( atoms[i][5].strip() ) + float( atoms[j][5].strip() ) ) / 2.}"
                        #epsilon (geometric awg of the two epsilons)
                        tmp_nonbond_params = tmp_nonbond_params + " " + f"{( float( atoms[i][6].strip() ) * float( atoms[j][6].strip() ) ) ** 0.5}"
                        #newline
                        tmp_nonbond_params = tmp_nonbond_params + "\n"

                        #append to the cumulative ones
                        nonbond_params.append(tmp_nonbond_params)

                nonbond_params.append("\n\n")

                #add the nonbond_params in untouched_itp 
                untouched_itp[i:i] = nonbond_params

                break


            if is_atomtype and (untouched_itp[i].strip() != "") and (untouched_itp[i].strip()[0] != ";"):

                line = untouched_itp[i].strip().split()

                tmp_DUM = ["DUM_" + line[0].strip(), line[1].strip(), line[2].strip(), "0.0", line[4].strip(), "0.0", "0.0", "\n"]
                
                tmp_DUM = " ".join(tmp_DUM)

                #every atom will be daclared with its DUM_ version
                untouched_itp[i] = untouched_itp[i].strip() + "\n" + tmp_DUM

                #keeping the normal atoms in order to do the nonbond_params later
                atoms.append(line)

        #itp for ligand gets created (only_ligand)
        ligand_creation_itp = untouched_itp

        
        is_atoms = False
        for i in range(len(ligand_creation_itp)):

            if ligand_creation_itp[i].strip().replace(" ", "") == "[atoms]" and is_atoms == False:

                is_atoms = True

            elif is_atoms and ligand_creation_itp[i].strip()[0] == "[":

                break

            if is_atoms and (ligand_creation_itp[i].strip() != "") and (ligand_creation_itp[i].strip()[0] != ";"):

                tmp_line = ligand_creation_itp[i].strip().split()

                #add DUM_ state a
                ligand_creation_itp[i] = tmp_line[0]

                ligand_creation_itp[i] = ligand_creation_itp[i] + " " + "DUM_" + tmp_line[1]

                ligand_creation_itp[i] = ligand_creation_itp[i] + " " + tmp_line[2]

                ligand_creation_itp[i] = ligand_creation_itp[i] + " " + tmp_line[3]

                ligand_creation_itp[i] = ligand_creation_itp[i] + " " + tmp_line[4]

                ligand_creation_itp[i] = ligand_creation_itp[i] + " " + tmp_line[5]

                ligand_creation_itp[i] = ligand_creation_itp[i] + " " + "0.0"

                ligand_creation_itp[i] = ligand_creation_itp[i] + " " + tmp_line[7]

                ligand_creation_itp[i] = ligand_creation_itp[i] + " " + tmp_line[2]

                ligand_creation_itp[i] = ligand_creation_itp[i] + " " + tmp_line[6]

                ligand_creation_itp[i] = ligand_creation_itp[i] + " " + tmp_line[7]

                ligand_creation_itp[i] = ligand_creation_itp[i] + "\n"

        #writing it in the right place
        write_on_files.write_file(lines = ligand_creation_itp, file_name = itp_filename)


    def _merge_ligand_solvent_files(self, fsdam_dir, solvent_gro):
        
        only_solvent = read_file.read_file(file_name = solvent_gro)

        #removing possible empty lines from end of the list
        lines_to_remove = []
        for i in range(len(only_solvent) - 1, 0, -1):

            if only_solvent[i].strip() == "":
                lines_to_remove.append(i)

            else:
                break

        if lines_to_remove != []:

            del only_solvent[lines_to_remove[-1]:lines_to_remove[0] + 1]


        solvent_atoms = int(only_solvent[1].strip())

        box = only_solvent[-1]

        #remove the useless lines
        del only_solvent[0:2]

        del only_solvent[-1]

        #add the solvent box to any gro file in fsdam_dir
        for item in os.listdir(fsdam_dir):

            if item[-4:] == ".gro":

                ligand_gro = read_file.read_file(fsdam_dir + "/" + item)

                ligand_gro[1] = f" {solvent_atoms + int(ligand_gro[1].strip())} \n"

                #removing possible empty lines from end of the list
                lines_to_remove = []
                for i in range(len(ligand_gro) - 1, 0, -1):

                    if ligand_gro[i].strip() == "":
                        lines_to_remove.append(i)

                    else:
                        break

                if lines_to_remove != []:

                    del ligand_gro[lines_to_remove[-1]:lines_to_remove[0] + 1]

                #temporarely remove the box line
                del ligand_gro[-1]

                residue_number = int(ligand_gro[-1][0:5].strip()) + 1

                atom_number = int(ligand_gro[-1][15:20]) + 1

                old_residue_number = int(only_solvent[0][0:5].strip())

                for line in only_solvent:

                    if int(only_solvent[0][0:5].strip()) != old_residue_number:

                        residue_number = residue_number + 1

                        old_residue_number = int(only_solvent[0][0:5].strip())

                    ligand_gro.append("{:>5}{}{:>5}{}".format(residue_number, line[5:10], atom_number, line[20:]))

                    atom_number = atom_number + 1

                ligand_gro.append(box)

                write_on_files.write_file(lines = ligand_gro, file_name = fsdam_dir + "/" + item)



    def _create_restart_configs(self,fsdam_dir, only_solvent_gro):
        """
        Private

        Makes the starting files for fs-dam and writes them in the RESTART dir
        """

        i = 0
        while os.path.exists(f"BATTERY{i}"):

            subprocess.run(f"echo System | {self.gromacs_path} trjconv -f BATTERY{i}/scaled0/traj_comp.xtc -s BATTERY{i}/scaled0/topol.tpr -pbc mol -b 80 -o {fsdam_dir}/BATTERY{i}_.gro -sep yes -ur compact", shell = True)
            
            i = i + 1

        self._merge_ligand_solvent_files(fsdam_dir = fsdam_dir, solvent_gro = only_solvent_gro)


    def _edit_top(self, ligand_top, only_solvent_top):

        only_solvent_lines = read_file.read_file(file_name = only_solvent_top)

        for i in range(len(only_solvent_lines) - 1, -1, -1):

            if only_solvent_lines[i].strip() != "":

                solvent_line = only_solvent_lines

                break

        ligand_lines = read_file.read_file(file_name = ligand_top)

        ligand_lines.append(solvent_line)

        write_on_files.write_file(lines = ligand_lines, file_name = ligand_top)

    
    def _prepare_input(self, ligand_resname, fsdam_dir, top_file):
        
        prepare_input = PrepareInputOnlyLigand(fsdam_dir = fsdam_dir,
                            ligand_resname = ligand_resname,
                            top_file = top_file,
                            cpus_per_node = self.cpus_per_node,
                            use_GPU = self.use_GPU,
                            GPU_per_node = self.GPU_per_node)

        prepare_input.execute()

    
    def execute(self):
        
        #move to the HREM directory
        if self.HREM_dir != os.getcwd():
            
            changed_dir = True

            old_dir = os.getcwd()

            os.chdir(self.HREM_dir)

        else:

            changed_dir = False

            old_dir = None


        fsdam_dir = "RESTART"

        os.makedirs(fsdam_dir, exist_ok=True)

        #lazily copy all the itp files in the new dir
        for i in os.listdir(os.getcwd()):

            if i[-4:] == ".itp":

                shutil.copy(i, fsdam_dir)

        #making some defaults
        useful_info = {
            f"ligand_resname" : "LIG",
            f"ligand_itp" : "LIG.itp",
            f"ligand_top" : "LIG.top",
            f"only_solvent_gro" : None,
            f"only_solvent_top" : None
        }

        lines = read_file.read_file(file_name = "important_info.dat")
        for line in lines:

            line = line.strip()
            line = line.split("=")
            useful_info[line[0].strip()] = line[1].strip()

        self._create_restart_configs(fsdam_dir = fsdam_dir, only_solvent_gro = useful_info["only_solvent_gro"])


        #creates the two new itp files (ligand creation and annihilation)
        self._create_itp_files(input_itp_filename = useful_info['ligand_itp'],
                            itp_filename = f"{fsdam_dir}/{useful_info['ligand_itp']}")

        
        #edit top file
        self._edit_top(ligand_top = useful_info["top_file"], only_solvent_top = useful_info["only_solvent_top"])

        #copy the protein ligand top in the protein_ligand dir
        shutil.copy(useful_info["top_file"], fsdam_dir)

        #prepares the input files for the 
        self._prepare_input(ligand_resname = useful_info["ligand_resname"],
                            fsdam_dir = fsdam_dir,
                            top_file = useful_info["top_file"])


        #go back to the old working directory (attention, if something fails the working dir will remain in the
        # HREM_dir!! (use try: except: in case needed))

        if changed_dir:

            os.chdir(old_dir)









class PrepareInputSuperClass(object):

    """
    It is an ancillary class for the fsdam class
    this is it's superclass
    """

    def __init__(self,
                fsdam_dir,
                ligand_resname,
                top_file,
                cpus_per_node = 64,
                use_GPU = False,
                GPU_per_node = 1):

        self.fsdam_dir = fsdam_dir

        self.ligand_resname = ligand_resname

        self.top_file = top_file

        self.cpus_per_node = cpus_per_node

        self.use_GPU = use_GPU

        self.GPU_per_node = GPU_per_node

    def _make_transitionQ_mdp(self, mdp_file_name):

        mdp_lines = []

        write_on_files.write_file(lines = ["\n".join(mdp_lines)], file_name = mdp_file_name)


    def _make_transitionVdW_mdp(self, mdp_file_name):

        mdp_lines = []

        write_on_files.write_file(lines = ["\n".join(mdp_lines)], file_name = mdp_file_name)


    def _write_make_tpr_script(self, number_of_dirs, mdp_file, gro_file, tpr_file, script_filename):

        script = [
            "#!/bin/bash\n",
            f"for i in {{1..{number_of_dirs}}}; do\n",
            f"    gmx grompp -f {mdp_file} -c fsdam_{number_of_dirs}/{gro_file} -p {self.top_file} -maxwarn 100 -o fsdam_{number_of_dirs}/{tpr_file}\n",
            "done\n"
        ]

        write_on_files.write_file(lines = script, file_name = script_filename)


    def _make_workload_manager_input(self, number_of_dirs, tpr_file):

        prefix = tpr_file[:-4]

        multidir = ""
        for i in range(number_of_dirs):

            multidir = multidir + f" fsdam_{i}"

        mpirun_string = "\nmpirun -np {} gmx_mpi mdrun -deffnm {} -s {} -multidir {} "

        if self.use_GPU:
            #one GPU per fsdam

            nodes = math.ceil(number_of_dirs / self.GPU_per_node)

            tasks = number_of_dirs

            tasks_per_node = self.GPU_per_node

            cpus_per_task = math.ceil(self.cpus_per_node / tasks_per_node)

            GPUs = tasks

            mpirun_string = mpirun_string.format(cpus_per_task * tasks, prefix, tpr_file, multidir)

        else:
            #8 cpu per fsdam

            tasks = number_of_dirs

            cpus_per_task = 8

            tasks_per_node = math.floor(self.cpus_per_node / cpus_per_task)

            nodes = math.ceil(number_of_dirs / tasks_per_node)

            GPUs = None

            mpirun_string = mpirun_string.format(cpus_per_task * tasks, prefix, tpr_file, multidir)


        #slurm input
        slurm = workload_managers.SlurmHeader(nodes = nodes,
                                            tasks = tasks,
                                            tasks_per_node = tasks_per_node,
                                            cpus_per_task = cpus_per_task,
                                            GPUs = GPUs,
                                            wall_time = "24:00:00",
                                            output = f"{prefix}_stdout.out",
                                            error = f"{prefix}_stderr.err",
                                            account_name = None,
                                            partition_name = None)

        slurm_header = slurm.execute()

        slurm_header.append(mpirun_string)

        slurm_header = ["\n".join(slurm_header)]

        write_on_files.write_file(lines = slurm_header, file_name = self.fsdam_dir + "/" + f"{prefix}_input.slr")


        #pbs input
        pbs = workload_managers.PBSHeader(nodes = nodes,
                                            tasks = tasks,
                                            tasks_per_node = tasks_per_node,
                                            cpus_per_task = cpus_per_task,
                                            GPUs = GPUs,
                                            wall_time = "24:00:00",
                                            output = f"{prefix}_stdout.out",
                                            error = f"{prefix}_stderr.err",
                                            account_name = None,
                                            partition_name = None)

        pbs_header = pbs.execute()

        pbs_header.append(mpirun_string)

        pbs_header = ["\n".join(slurm_header)]

        write_on_files.write_file(lines = slurm_header, file_name = self.fsdam_dir + "/" + f"{prefix}_input.pbs")



    def execute(self):

        #creates the mdp files
        self._make_transitionQ_mdp(self.fsdam_dir + "/" + "transitionQ.mdp")
        self._make_transitionVdW_mdp(self.fsdam_dir + "/" + "transitionVdW.mdp")

        #creates 200 directories for the fsdam, and moves 200 gro files there
        #choses randomly 200 gro files out of the available ones
        from random import shuffle
        j = 0
        for i in shuffle(os.listdir(self.fsdam_dir)):

            if i[-4:] == ".gro":

                os.makedirs(self.fsdam_dir + "/" + f"fsdam_{j}", exist_ok = True)

                shutil.move(self.fsdam_dir + "/" + i, self.fsdam_dir + "/" + f"fsdam_{j}" + "/" + "start_fsdam.gro")

                if j > 200:
                    break

                j = j + 1

                
        #creates the bash scripts to create the needed tpr files
        self._write_make_tpr_script(number_of_dirs = j,
                                mdp_file = "transitionQ.mdp",
                                gro_file = "start_fsdam.gro",
                                tpr_file = "transitionQ.tpr",
                                script_filename = self.fsdam_dir + "/" + "MAKE_TPR_transitionQ.sh")

        self._write_make_tpr_script(number_of_dirs = j,
                                mdp_file = "transitionVdW.mdp",
                                gro_file = "transitionQ.gro",
                                tpr_file = "transitionVdW.tpr",
                                script_filename = self.fsdam_dir + "/" + "MAKE_TPR_transitionVdW.sh")

        #make the inputs for both Q and VdW
        self._make_workload_manager_input(number_of_dirs = j, tpr_file = "transitionQ.tpr")

        self._make_workload_manager_input(number_of_dirs = j, tpr_file = "transitionVdW.tpr")



class PrepareInputProteinLigand(PrepareInputSuperClass):

    def _make_transitionQ_mdp(self, mdp_file_name):

        mdp_lines = [
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
            "dt                       = 0.001",
            "nsteps                   = 1000000",
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
            "nstxout                  = 10000",
            "nstvout                  = 10000",
            "nstfout                  = 10000",
            "; Output frequency for energies to log file and energy file",
            "nstlog                   = 1000",
            "nstcalcenergy            = 100",
            "nstenergy                = 1000",
            "; Output frequency and precision for .xtc file",
            "nstxtcout                = 2000",
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
            "morse                    = no",
            "",
            "; Free energy control stuff",
            "free-energy              = yes",
            "init-lambda              = 1",
            "delta-lambda             = -1e-6",

            f"couple-moltype           = {self.ligand_resname}",

            "couple-lambda0           =vdw",
            "couple-lambda1           =vdw-q",
            "sc-alpha                 = 0.3",
            "sc-coul                  = yes",
            "sc-sigma                 = 0.25",
            "sc-power                 = 1"
        ]

        write_on_files.write_file(lines = ["\n".join(mdp_lines)], file_name = mdp_file_name)


    def _make_transitionVdW_mdp(self, mdp_file_name):

        mdp_lines = [
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
            "dt                       = 0.001",
            "nsteps                   = 1000000",
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
            "nstxout                  = 10000",
            "nstvout                  = 10000",
            "nstfout                  = 10000",
            "; Output frequency for energies to log file and energy file",
            "nstlog                   = 1000",
            "nstcalcenergy            = 100",
            "nstenergy                = 1000",
            "; Output frequency and precision for .xtc file",
            "nstxtcout                = 2000",
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
            "morse                    = no",
            "",
            "; Free energy control stuff",
            "free-energy              = yes",
            "init-lambda              = 1",
            "delta-lambda             = -1e-6",

            f"couple-moltype           = {self.ligand_resname}",

            "couple-lambda0           =none",
            "couple-lambda1           =vdw",
            "sc-alpha                 = 0.3",
            "sc-coul                  = yes",
            "sc-sigma                 = 0.25",
            "sc-power                 = 1"
        ]

        write_on_files.write_file(lines = ["\n".join(mdp_lines)], file_name = mdp_file_name)






class PrepareInputOnlyLigand(PrepareInputSuperClass):

    def _make_transitionQ_mdp(self, mdp_file_name):

        #TO DO
        mdp_lines = []

        write_on_files.write_file(lines = ["\n".join(mdp_lines)], file_name = mdp_file_name)


    def _make_transitionVdW_mdp(self, mdp_file_name):

        #TO DO
        mdp_lines = []

        write_on_files.write_file(lines = ["\n".join(mdp_lines)], file_name = mdp_file_name)