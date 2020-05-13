######################################################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock                                        #
#                                                                                    #
# This software is open-source and is distributed under the                          #
# GNU Affero General Public License v3 (agpl v3 license)                             #
#                                                                                    #
# A copy of the license must be included with any copy of the program or part of it  #
######################################################################################

"""
This file contains the classes and functions needed to get the input
"""
import importlib_resources

from HPC_Drug.auxiliary_functions import path as program_path

class GetInput(object):
    """Generic class to get input"""

    def __init__(self, *args, **kwargs):
        self.input = args
        self.dictionary = kwargs

class GetFile(GetInput):
    """Class to get filename stored
    any extra argument will be stored in self.input or self.dictionary"""

    def __init__(self, filename = None, *args, **kwargs):
        self.filename = filename
        if type(self.filename) != str:
            raise TypeError(f"filename must be a string, not a {type(self.filename)} type")
        super().__init__(args, kwargs)


class ParseInputFromFile(GetFile):
    """Class to get the input file to start the program
    Will consider only expression in form key = value , the rest is considered comment
    Not kase-sensitive
    It returns a dictionary
    
    :param:filename
    the name of the input file or it's absolute path"""

    

    def __init__(self, filename):
        super().__init__(filename)
        
        self.possible_keys = ('protein',
                            'ligand',
                            'protein_filetype',
                            'ligand_elaboration_program',
                            'ligand_elaboration_program_path',
                            'local',
                            'filepath',
                            'Protein_model',
                            'Protein_chain',
                            'ligand_in_protein',
                            'ph',
                            'repairing_method',
                            'MD_program',
                            'MD_program_path',
                            'protein_prm_file',
                            'protein_tpg_file',
                            'solvent_pdb',
                            'kind_of_processor',
                            'number_of_cores_per_node',
                            'residue_substitution',
                            'use_gpu',
                            'gpu_per_node')

        
        self.input_variables = self.read_input()

        
    def _create_input_dict(self):
        input_dict = {}

        for key in self.possible_keys:
            input_dict[key] = None
        
        return input_dict

    def read_input(self):
        """It returns a dictionary
    
        :param:filename
        the name of the input file or it's absolute path"""

        input_variables = self._create_input_dict()

        with open(self.filename) as f:

            for line in f:

                line = line.strip()

                # the '#' identifies a comment
                if len(line) == 0:
                    continue
                elif line[0] == '#':
                    continue

                line = line.split('=')

                if len(line) == 2:

                    key, value = line

                    key = key.strip()
                    value = value.strip()

                    if key in self.possible_keys:
                        input_variables[key] = value

                    else:
                        raise ValueError('InvalidInputKey ', key)
                
                else:
                    raise ValueError(f"Input must be key = value not like this: '{line}'")

        input_variables = self._refine_input(input_variables)    
        
        return input_variables

    def _refine_input(self, input_variables = None):

        """Makes the needed casts from string, defines some None to default ecc"""

        # Protein_model must be an integer
        if input_variables['Protein_model'] != None:
            input_variables['Protein_model'] = int(input_variables['Protein_model'])
        #default is zero
        elif input_variables['Protein_model'] == None:
            input_variables['Protein_model'] = 0

        #The standar chain is 'A'
        if input_variables['Protein_chain'] == None:
            input_variables['Protein_chain'] = 'A'
        #need it uppercase
        else:
            input_variables['Protein_chain'] = input_variables['Protein_chain'].upper().strip()

        #ph must be a float
        if input_variables['ph'] != None:
            input_variables['ph'] = float(input_variables['ph'])
        
        #sets primadorac as the default ligand elaboration program and tries to guess where the executable is
        if input_variables['ligand_elaboration_program'] == None:
            input_variables['ligand_elaboration_program'] = 'primadorac'
        
        if input_variables['ligand_elaboration_program_path'] == None:
            input_variables['ligand_elaboration_program_path'] = program_path.absolute_programpath(program = '~/ORAC/trunk/tools/primadorac/primadorac.bash')
        else:
            input_variables['ligand_elaboration_program_path'] = program_path.absolute_programpath(program = input_variables['ligand_elaboration_program_path'])

        #absolute path of the program path
        input_variables['MD_program_path'] = program_path.absolute_programpath(program = input_variables['MD_program_path'])


        if input_variables['MD_program'] == None:
            input_variables['MD_program'] = 'gromacs'

        if input_variables['MD_program'] == 'orac':
            if input_variables['protein_tpg_file'] == None:
                with importlib_resources.path('HPC_Drug.lib', 'amber99sb-ildn.tpg') as path:
                    input_variables['protein_tpg_file'] = str(path.resolve())
            
            if input_variables['protein_prm_file'] == None:
                with importlib_resources.path('HPC_Drug.lib', 'amber99sb-ildn.prm') as path:
                    input_variables['protein_prm_file'] = str(path.resolve())

            if input_variables['solvent_pdb'] == None:
                with importlib_resources.path('HPC_Drug.lib', 'water.pdb') as path:
                    input_variables['solvent_pdb'] = str(path.resolve())
        
        elif input_variables['MD_program'] == 'gromacs':
            if input_variables['protein_tpg_file'] == None:
                input_variables['protein_tpg_file'] = 'amber99sb-ildn'

            if input_variables['solvent_pdb'] == None:
                input_variables['solvent_pdb'] = 'spce'


        #set a default for processor (skylake)
        if input_variables['kind_of_processor'] == None:
            input_variables['kind_of_processor'] = 'skylake'
        
        else:
            #in this way I don't have to worry about strange input formats
            input_variables['kind_of_processor'] = input_variables['kind_of_processor'].lower().strip()

        #the number of cores per node is usually 64
        if input_variables['number_of_cores_per_node'] == None:
            input_variables['number_of_cores_per_node'] = 64
        
        else:
            #cast string to integer
            input_variables['number_of_cores_per_node'] = int(input_variables['number_of_cores_per_node'].strip())

        
        if input_variables['residue_substitution'] == None:
            input_variables['residue_substitution'] = 'standard'

        if input_variables['use_gpu'] == None:
            input_variables['use_gpu'] = 'auto'

        if input_variables['gpu_per_node'] is None:
            input_variables['gpu_per_node'] = 1
        
        else:
            input_variables['gpu_per_node'] = int(input_variables['gpu_per_node'])
        
        return input_variables

        



        
