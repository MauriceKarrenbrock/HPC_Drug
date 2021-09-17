######################################################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains some high level classes to rename the residues acording to the FF and
the ff options chosen

And also the concrete classes that implement it
"""

import warnings
from HPC_Drug import important_lists

class ResidueRenamer(object):
    """
    This class choses the right classes to use
    acording to the input (factory)
    """

    def __init__(self,
                Protein,
                forcefield='amber',
                substitution = "standard",
                ph = 7.0):

        if forcefield == 'amber':
            
            if substitution == "standard":

                self.residue_substitutor = AmberResidueRenamerStandardSubstitution(Protein = Protein, ph = ph)

            elif substitution == "custom_zinc":

                self.residue_substitutor = AmberResidueRenamerCustomZincSubstitution(Protein = Protein, ph = ph)

            else:
                raise NotImplementedError(f"{substitution} is not an implemented substitution method (remember that it is case sensitive)")


        else:
            raise NotImplementedError(f"The forcefield {forcefield} is not an implemented option, remember that the input is case sensitive")


        self.Protein = Protein

    def execute(self):

        self.Protein = self.residue_substitutor.execute()

        return self.Protein


        

class SpecificResidueRenamerSuperclass(object):
    """
    This is the superclass for all the specific
    residue renaming classes (forcefield specific, substitution type specific etc)
    """

    def __init__(self):

        self.Protein = None

        raise Exception("SpecificResidueRenamerSuperclass is a superclass, you cannot instantiate it")
    
    def _check_if_residue_binds_metal_with_two_atoms_instead_of_one(self, residue, metal = "ZN", tollerance = 0.01):

        import copy

        structure = copy.deepcopy(self.Protein.structure)

        for model in structure:
            for chain in model:
                for other_residue in chain:
                    for other_atom in other_residue:

                        if other_atom.name.upper() == metal:

                            shortest_distance = residue[self.Protein.substitutions_dict[residue.id[1]][1].strip().upper()] - other_atom

                            if shortest_distance < 3.:

                                for atom in residue:
                                    if atom.name.upper() != self.Protein.substitutions_dict[residue.id[1]][1].strip().upper():

                                        tmp_distance = atom - other_atom

                                        if tmp_distance - shortest_distance <= tollerance:

                                            #there are 2 atoms in the residue that bind the metallic ion
                                            return True

        #only one atom in the residue binds the metallic ion
        return False



    def _hook(self, residue):
        """
        PRIVATE
        it is an empty method that will
        be implemented by the subclasses
        """

        raise NotImplementedError("This hook method has not been impemented")


    def _iterate_residues(self, structure):
        """
        PRIVATE
        Iterates through a given biopython structure residue by residue
        and executes a _hook method (self._hook(residue)) that will
        be implemented by subclasses

        returns the new structure
        """

        for model in structure:

            for chain in model:

                for residue in chain:

                    residue = self._hook(residue)

        return structure

    def execute(self):

        raise NotImplementedError("This method has not been impemented")


#######################################
# AMBER
########################################

class AmberResidueRenamerStandardSubstitution(SpecificResidueRenamerSuperclass):

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

        return residue


    def execute(self):

        self.Protein.update_structure(struct_type = "biopython")

        self._iterate_residues(structure = self.Protein.structure)

        self.Protein.write(struct_type = 'biopython')

        return self.Protein



class AmberResidueRenamerCustomZincSubstitution(SpecificResidueRenamerSuperclass):

    def __init__(self, Protein, ph = 7.0):

        self.Protein = Protein

        self.ph = ph

        self.standard_substitutor = AmberResidueRenamerStandardSubstitution(Protein = self.Protein, ph = self.ph)


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

        return residue


    def execute(self):

        self.Protein = self.standard_substitutor.execute()

        self.Protein.update_structure(struct_type = "biopython")

        self._iterate_residues(structure = self.Protein.structure)

        self.Protein.write(struct_type = 'biopython')

        return self.Protein

########################################
# END AMBER
#########################################
