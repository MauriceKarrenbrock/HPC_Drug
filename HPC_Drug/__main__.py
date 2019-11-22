# This is the main file
import sys
import get_input
import pipelines
import os



try:
    input_file_name = sys.argv[1]
except IndexError as err:
    raise IndexError("Needs a filename as input")

try:
    input_variables = get_input.ParseInputFromFile(input_file_name).input_variables
except FileNotFoundError as err:
    raise FileNotFoundError("Did not find the file")
except ValueError as err:
    raise ValueError(err.args)

pipeline = pipelines.choose_pipeline(**input_variables)

pipeline.execute()