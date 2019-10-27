# This is the main file
import sys
import get_input
import structures
import pipelines
import os

try:
    input_file_name = sys.argv[1]
except IndexError as err:
    raise IndexError("Needs a filename as input")

try:
    input_variables = get_input.GetInputFromFile(input_file_name).input_variables
except FileNotFoundError as err:
    raise FileNotFoundError("Did not find the file")
except ValueError as err:
    raise ValueError(err.args)

# Downloads the protein structure file, if not already present
structures.Protein.download_protein_structure(protein_id = input_variables['protein'],\
    file_type = input_variables['protein_filetype'], pdir = os.getcwd())

Pipeline = ChoosePipeline(input_variables['protein'], input_variables['protein_filetype'], input_variables['ligand'])