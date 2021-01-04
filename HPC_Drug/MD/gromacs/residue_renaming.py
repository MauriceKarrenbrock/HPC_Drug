######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

from HPC_Drug import important_lists
from HPC_Drug.MD import residue_renaming

import warnings


class GromacsResidueRenamerStandardSubstitution(residue_renaming.SpecificResidueRenamerSuperclass):

    def __init__(self, Protein, ph = 7.0):

        self.Protein = Protein

        self.ph = ph

    def _rename_hist(self, residue):

        res_id = residue.id[1]

        if res_id in self.Protein.substitutions_dict.keys():

            if self.Protein.substitutions_dict[res_id][2] not in important_lists.metals:

                warnings.warn(f"{self.Protein.substitutions_dict[res_id][2]} is not a known metal, thus residue nr {res_id} could have been renamed wrongly")

            if self.Protein.substitutions_dict[res_id][1] == 'NE2':
                residue.resname = 'HID'
        
            elif self.Protein.substitutions_dict[res_id][1] == 'ND1':
                residue.resname = 'HIE'
            
            else:
                raise NotImplementedError(f"{self.Protein.substitutions_dict[res_id][1]} is not an implemented atom name for histidine")

        else:

            if self.ph < 6.0:
                residue.resname = 'HIS'
            
            else:
                residue.resname = 'HID'

        return residue

                    
    def _rename_cyst(self, residue):

        res_id = residue.id[1]

        if res_id in self.Protein.substitutions_dict.keys():

            if self.Protein.substitutions_dict[res_id][2].lower().strip() == "disulf":

                residue.resname = 'CYX'

            else:

                if self.Protein.substitutions_dict[res_id][2] not in important_lists.metals:

                    warnings.warn(f"{self.Protein.substitutions_dict[res_id][2]} is not a known metal, thus residue nr {res_id} could have been renamed wrongly")


                residue.resname = "CYM"


        return residue


    def _hook(self, residue):

        resname = residue.resname.upper()

        if resname in important_lists.hist_resnames:
            
            residue = self._rename_hist(residue = residue)

        elif resname in important_lists.cyst_resnames:

            residue = self._rename_cyst(residue = residue)


        elif resname in ("GLU", "ASP"):

            #no renaming needed in amber ff
            pass

        else:
            warnings.warn(f"{resname} renaming is not implemented, going on without changing it")

        return residue


    def execute(self):

        self.Protein.update_structure(struct_type = "biopython")

        self._iterate_residues(structure = self.Protein.structure)

        self.Protein.write(struct_type = 'biopython')

        return self.Protein



class GromacsResidueRenamerCustomZincSubstitution(residue_renaming.SpecificResidueRenamerSuperclass):

    def __init__(self, Protein, ph = 7.0):

        self.Protein = Protein

        self.ph = ph

        self.standard_substitutor = GromacsResidueRenamerStandardSubstitution(Protein = self.Protein, ph = self.ph)


    def _rename_glu(self, residue):

        residue.resname = "GLZ"

        if not self._check_if_residue_binds_metal_with_two_atoms_instead_of_one(residue = residue, metal = "ZN"):

            for atom in residue:

                if atom.name == self.Protein.substitutions_dict[residue.id[1]][1].strip().upper():

                    atom.fullname = " OZ "

                    break

        return residue


    def _rename_asp(self, residue):

        residue.resname = "ASZ"

        if not self._check_if_residue_binds_metal_with_two_atoms_instead_of_one(residue = residue, metal = "ZN"):

            for atom in residue:

                if atom.name == self.Protein.substitutions_dict[residue.id[1]][1].strip().upper():

                    atom.fullname = " OZ "

                    break

        return residue


    def _hook(self, residue):

        res_id = residue.id[1]

        if res_id in self.Protein.substitutions_dict.keys():
            
            if self.Protein.substitutions_dict[res_id][2].strip().upper() == "ZN":

                resname = residue.resname.upper()

                if resname == "HID":
                    
                    residue.resname = "HDZ"
                    
                elif resname == "HIE":
                    
                    residue.resname = "HEZ"

                elif resname == "CYM":

                    residue.resname = "CYZ"


                elif resname == "GLU":

                    residue = self._rename_glu(residue = residue)

                elif resname == "ASP":

                    residue = self._rename_asp(residue = residue)

                else:
                    warnings.warn(f"{resname} renaming is not implemented, going on without changing it")

        return residue


    def execute(self):

        self.Protein = self.standard_substitutor.execute()

        self.Protein.update_structure(struct_type = "biopython")

        self._iterate_residues(structure = self.Protein.structure)

        self.Protein.write(struct_type = 'biopython')

        return self.Protein


