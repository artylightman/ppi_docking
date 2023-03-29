# author : artemis lightman
# date created : jan 2, 2023
# last modified: jan 2, 2023

# command line arguments
    # none

# dependencies:
    # arty.py

# input: pdb decoys in txt form
# output: 
    # models.txt - list of model file names

##############################################################################
# creates list of model file names
##############################################################################

import glob
from datetime import datetime
import argparse
import arty

dt_start = datetime.now()

parser = argparse.ArgumentParser(description = "this script creates a list of model file names")

args = parser.parse_args()

filenames = glob.glob("*_complex.txt")

output = open("models.txt", "w")

for item in filenames:

	output.write(item[:-12] + "\n")

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y")

str_runtime = arty.get_runtime(dt_start, dt)

print ("print_model_names.py finished at " + str_dt + " (runtime " + str_runtime + ")")