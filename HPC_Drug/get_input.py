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

    possible_keys = ('protein', 'ligand',\
        'protein_filetype', 'ligand_elaboration_program',\
            'local', 'filepath', 'PDB_model')

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

                line = line.replace(" ", "")
                line = line.strip()

                # the '#' identifies a comment
                if len(line) == 0:
                    continue
                elif line[0] == '#':
                    continue

                line = line.split('=')

                if len(line) == 2:

                    key, value = line

                    if key in self.possible_keys:
                        input_variables[key] = value

                    else:
                        raise ValueError('InvalidInputKey ', key)
                
                else:
                    raise ValueError(f'Imput must be key = value not like this: "{line}"')

            # PDB_model must be an integer
            if input_variables['PDB_model'] != None:
                input_variables['PDB_model'] = int(input_variables['PDB_model'])

            return input_variables


        
