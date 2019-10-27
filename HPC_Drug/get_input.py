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


class GetInputFromFile(GetFile):
    """Class to get the input file to start the program
    Will consider only expression in form key = value , the rest is considered comment
    Not kase-sensitive
    
    It returns a dictionary """

    possible_keys = ('protein', 'ligand', 'protein_filetype', 'ligand_elaboration_program')

    def __init__(self, filename):
        super().__init__(filename)
        self.input_variables = self.read_input()



    def create_input_dict(self):
        input_dict = {}

        for key in self.possible_keys:
            input_dict[key] = None
        
        return input_dict

    def read_input(self):

        input_variables = self.create_input_dict()

        with open(self.filename) as f:
            lines = f.readlines()

            for line in lines:
                line = line.split('=')

                if len(line) == 2:

                    key, value = line

                    key = key.strip()
                    value = value.strip()
                    key = key.lower()
                    value = value.lower()

                    if key in self.possible_keys:
                        input_variables[key] = value

                    else:
                        raise ValueError('InvalidInputKey ', key)
                
            return input_variables


        
