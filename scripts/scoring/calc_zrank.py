# author : artemis lightman
# date created : jan 2, 2023
# last modified: jan 2, 2023

# command line arguments
	# t/target - indicates whether structures are targets
	# n - number of structures
	# id - pdb id of target
	# chains - chains in target
	# em - indicates that structures are energy minimized (affects filename)
	# r - indicates that structures are rosetta relaxed (affects filename)

# dependencies:
	# arty.py, zrank

# input: none
# output: zrank score of structures

##############################################################################
# calculates zrank score of structures
##############################################################################

import os
from datetime import datetime
import argparse
import arty
import time

dt_start = datetime.now()

parser = argparse.ArgumentParser(description = "this script calculates zrank scores")

parser.add_argument("-t", "--target", help = "indicates that the structures are targets", action = "store_true")
parser.add_argument("-n", help = "number of structures", default = 500, type = int)
parser.add_argument("-id", help = "pdb id of target", default = "1acb")
parser.add_argument("-em", help = "indicates that the structures are energy minimized (affects filename)", action = "store_true")
parser.add_argument("-r", help = "indicates that the structures are rosetta relaxed (affects filename)", action = "store_true")

args = parser.parse_args()

pdb_names = []
chains = []

if args.target:

	with open("pdb_chains.txt") as input:

		for line in input:

			arr = line.split()

			pdb_id = arr[0]

			chains_target = arr[1]

			chains.append(chains_target)

			if args.em:

				filename = "em_" + pdb_id

			elif args.r:

				filename = "relaxed_" + pdb_id

			else:

				filename = pdb_id

			pdb_names.append(filename)

else:

	with open("../pdb_chains.txt") as input:

		for line in input:

			if line.startswith(args.id):

				chains_model = line.split()[1]

	for i in range(1, int(args.n) + 1):

		if args.em:

			filename = "em_" + args.id + "_model_" + str(i)

		elif args.r:

			filename = "relaxed_" + args.id + "_model_" + str(i)

		else:

			filename = args.id + "_model_" + str(i)

		print (filename)

		if os.path.exists(filename + "_complex_H.pdb"):

			chains.append(chains_model)
			pdb_names.append(filename)

for index, structure in enumerate(pdb_names):

	target = structure

	if args.em:

		chain_one = "A"
		chain_two = "B"

	else:

		chain_one = chains[index][0]
		chain_two = chains[index][1]

	######## add TER lines for zrank to run ########

	ter_added = open(target + "_ter.pdb", "w")
	chain_one_file = target + "_" + chain_one + "_H.pdb"
	chain_two_file = target + "_" + chain_two + "_H.pdb"

	with open(chain_one_file) as input:

		for line in input:

			if line.startswith("ATOM"):

				ter_added.write(line)

	ter_added.write("TER" + "\n")

	with open(chain_two_file) as input:

		for line in input:

			if line.startswith("ATOM"):

				ter_added.write(line)

zrank_list = open("zrank_list.txt", "w")

for n in pdb_names:

	filename = n + "_ter.pdb"

	if os.path.exists(filename):

		zrank_list.write(filename + "\n")

zrank_list.close()

time.sleep(5) ## allows list to be written before zrank runs

os.system("./zrank zrank_list.txt")

os.system("rm *_ter.pdb")

if args.target:

	os.system("mv zrank_list.txt.zr.out zrank_scores_targets.txt")

else:

	os.system("mv zrank_list.txt.zr.out zrank_scores_" + args.id + ".txt")

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y")

str_runtime = arty.get_runtime(dt_start, dt)

print ("calc_zrank.py finished at " + str_dt + " (runtime " + str_runtime + ")")

