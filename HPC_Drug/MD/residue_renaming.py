######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains some high level classes to rename the residues acording to the MD program and
the ff options chosen
"""

class ResidueRenamer(object):
    """
    This class choses the right classes to use
    acording to the input (factory)
    """

    def __init__(self,
                Protein,
                MD_program,
                substitution = "standard",
                ph = 7.0):

        if MD_program == "orac":

            from HPC_Drug.MD.orac import residue_renaming
            
            if substitution == "standard":

                self.residue_substitutor = residue_renaming.OracResidueRenamerStandardSubstitution(Protein = Protein, ph = ph)

            elif substitution == "custom_zinc":

                self.residue_substitutor = residue_renaming.OracResidueRenamerCustomZincSubstitution(Protein = Protein, ph = ph)

            else:
                raise NotImplementedError(f"{substitution} is not an implemented substitution method (remember that it is case sensitive)")

        elif MD_program == "gromacs":

            from HPC_Drug.MD.gromacs import residue_renaming
            
            if substitution == "standard":

                self.residue_substitutor = residue_renaming.GromacsResidueRenamerStandardSubstitution(Protein = Protein, ph = ph)

            elif substitution == "custom_zinc":

                self.residue_substitutor = residue_renaming.GromacsResidueRenamerCustomZincSubstitution(Protein = Protein, ph = ph)

            else:
                raise NotImplementedError(f"{substitution} is not an implemented substitution method (remember that it is case sensitive)")


        else:
            raise NotImplementedError(f"MD_program = {MD_program} is not an implemented option, remember that the input is case sensitive")


        self.Protein = Protein

    def execute(self):

        self.Protein = self.residue_substitutor.execute()

        return self.Protein


        

class SpecificResidueRenamerSuperclass(object):
    """
    This is the superclass for all the specific
    residue renaming classes (MDprogram specific, substitution type specific etc)
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
