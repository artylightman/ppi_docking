# author : artemis lightman
# date created : jan 3, 2023
# last modified: jan 3, 2023

# command line arguments
	# n - number of decoys per target

# dependencies:
	# arty.py, generate_decoys.py, targets.txt

# input: file with list of target names
# output: decoys for all targets

##############################################################################
# generates decoys for all targets
##############################################################################

import os
from datetime import datetime
import argparse
import arty

dt_start = datetime.now()

parser = argparse.ArgumentParser(description = "this script generates decoys for all targets using rosetta")

parser.add_argument("-n", type = int, help = "number of decoys to generate", default = 500)

args = parser.parse_args()

output = open("tasklist.sh", "w")

with open("targets.txt") as input:

	for line in input:

		arr = line.split()

		pdb_id = arr[0]

		chains = arr[1]

		output.write("python generate_decoys.py -id " + pdb_id + " -c " + chains + " -n " + str(args.n) + " \n")

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y ")

str_runtime = arty.get_runtime(dt_start, dt)

print ("generate_all.py finished at " + str_dt + " (runtime " + str_runtime + ")")