# author : artemis lightman
# date created : jan 2, 2023
# last modified: jan 2, 2023

# command line arguments
	# em - indicates that structures are energy minimized
	# id - pdb id of target

# dependencies:
	# arty.py 

# input: all pdb decoys for target
# output: ialign score output

##############################################################################
# calculates ialign scores
##############################################################################

import os
from datetime import datetime
import argparse
import arty

dt_start = datetime.now()

parser = argparse.ArgumentParser(description = "this script calculates iAlign scores")

parser.add_argument("-em", help = "indicates that the structures are energy minimized (affects chain IDs)", action = "store_true")
parser.add_argument("-id",  help = "pdb id of target", default = "1acb")

args = parser.parse_args()

pdb_id = args.id

os.system("cp ../../../targets/" + pdb_id + "_complex_H.pdb " + pdb_id + "_complex_H.pdb")

os.system("cp -r ../../../structure_calculations/ialign ialign")

os.system("cp *_complex_H.pdb ialign/pdb")

with open("../pdb_chains.txt") as input:

	for line in input:

		if line.startswith(pdb_id):

			chains = line.split()[1]

os.chdir("ialign/pdb")

# generate list of structures to run ialign/is-score

if args.em:

	os.system("python generate_list.py AB")

else:

	os.system("python generate_list.py " + chains)

os.system("../bin/isscore.pl -w scratch -l models.lst " + pdb_id + "_complex_H.pdb " + chains + " > ../" + pdb_id + "_ialign_output.txt")

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y")

str_runtime = arty.get_runtime(dt_start, dt)

print ("calc_ialign.py finished at " + str_dt + " (runtime " + str_runtime + ")")