# This file contains the classes and functions needed to get the input

class GetInput(object):
    """Generic class to get input"""

    def __init__(self, *args, **kwargs):
        self.input = args
        self.dictionary = kwargs

class GetFile(GetInput):
    """Class to get filename stored
    any extra argument will be stored in self.input or self.dictionary"""

    def __init__(self, filename, *args, **kwargs):
        self.filename = filename
        super().__init__(args, kwargs)


class ParseInputFromFile(GetFile):
    """Class to get the input file to start the program
    Will consider only expression in form key = value , the rest is considered comment
    Not kase-sensitive
    It returns a dictionary
    
    :param:filename
    the name of the input file or it's absolute path"""

    possible_keys = ('protein', 'ligand',
                    'protein_filetype',
                    'ligand_elaboration_program',
                    'local', 'filepath', 'Protein_model',
                    'ligand_in_protein', 'ph',
                    'repairing_method')

    def __init__(self, filename):
        super().__init__(filename)
        self.input_variables = self.read_input()

    def create_input_dict(self):
        input_dict = {}

        for key in self.possible_keys:
            input_dict[key] = None
        
        return input_dict

    def read_input(self):
        """It returns a dictionary
    
        :param:filename
        the name of the input file or it's absolute path"""

        input_variables = self.create_input_dict()

        with open(self.filename) as f:
            lines = f.readlines()

            for line in lines:

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
                    raise ValueError(f'Imput must be key = value not like this: "{line}"')

        input_variables = self.refine_input(input_variables)    
        
        return input_variables

    def refine_input(self, input_variables = None):

        """Makes the needed casts from string"""

        if type(input_variables) != dict:
            raise Exception('Need a dictionary')

        # Protein_model must be an integer
        if input_variables['Protein_model'] != None:
            input_variables['Protein_model'] = int(input_variables['Protein_model'])

        #ph must be a float
        if input_variables['ph'] != None:
            input_variables['ph'] = float(input_variables['ph'])
        
        return input_variables

        



        
