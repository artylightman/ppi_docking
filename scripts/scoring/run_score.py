# author : artemis lightman
# date created : jan 2, 2023
# last modified: jan 2, 2023

# command line arguments
    # n - number of pdb structures
    # em - indicates structures are energy minimized
    # r - indicates structures are rosetta relaxed

# dependencies:
    # arty.py

# input: pdb structures and packing tasklist
# output: none

##############################################################################
# runs hdock, zrank, and rosetta scoring for all decoys
##############################################################################

from datetime import datetime
import argparse
import arty

dt_start = datetime.now()

parser = argparse.ArgumentParser(description = "this script runs hdock, zrank, and rosetta scoring for all decoys")
parser.add_argument("-n", help = "number of PDB structures", type = int, default = 500)
parser.add_argument("-em", help = "indicates that the structures are energy minimized (affects filename)", action = "store_true")
parser.add_argument("-r", help = "indicates that the structures are rosetta relaxed (affects filename)", action = "store_true")

args = parser.parse_args()

pdb_names = []
chains = []

with open("pdb_chains.txt") as input:

	for line in input:

		arr = line.split()
		pdb_name = arr[0]
		chains_target = arr[1]
		pdb_names.append(pdb_name)

output = open("tasklist_scoring.sh", "w")

for index, item in enumerate(pdb_names):

	if args.em:

		output.write("python run_score_target.py -n " + str(args.n) + " -id " + item + " -em" + "\n")

	elif args.r:

		output.write("python run_score_target.py -n " + str(args.n) + " -id " + item + " -r" + "\n")

	else:

		output.write("python run_score_target.py -n " + str(args.n) + " -id " + item + "\n")

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y")

str_runtime = arty.get_runtime(dt_start, dt)

print ("run_score.py finished at " + str_dt + " (runtime " + str_runtime + ")")


