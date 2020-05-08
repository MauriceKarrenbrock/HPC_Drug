######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

import os.path
import subprocess
import copy
import shutil

from HPC_Drug.auxiliary_functions import path
from HPC_Drug.files_IO import read_file, write_on_files
from HPC_Drug.MD.gromacs import gro2pdb
from HPC_Drug.PDB import biopython, prody
from HPC_Drug import orient


class FSDAMInputProteinLigand(object):

    def __init__(self,
                HREM_dir = os.getcwd(),
                delta_distance = 4.,
                gromacs_path = "gmx"):

        #the root directory of the HREM
        #default is working directory
        self.HREM_dir = HREM_dir
        
        self.delta_distance = delta_distance

        self.gromacs_path = path.absolute_programpath(program = gromacs_path)


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

        
    def _prepare_input(self):
        pass

    
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
        self._prepare_input()


        #go back to the old working directory (attention, if something fails the working dir will remain in the
        # HREM_dir!! (use try: except: in case needed))

        if changed_dir:

            os.chdir(old_dir)


        



class FSDAMInputOnlyLigand(object):

    def __init__(self):
        pass


    def _create_itp_files(self, input_itp_filename, ligand_annihilation_itp_filename, ligand_creation_itp_filename):


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
        write_on_files.write_file(lines = ligand_creation_itp, file_name = ligand_creation_itp_filename)
    
    def execute(self):
        pass

    