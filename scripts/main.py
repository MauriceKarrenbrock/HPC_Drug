# This is the main file
import sys
from HPC_Drug import get_input
from HPC_Drug import pipelines
import os

def main():
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

if __name__ == "__main__":

    main()