# author : artemis lightman
# date created : jan 2, 2023
# last modified: jan 2, 2023

# command line arguments
	# t/target - indicates whether structures are targets
	# n - number of structures
	# id - pdb id of target
	# em - indicates that the structures are energy minimized (affects filename)
	# r - indicates that the structures are rosetta relaxed (affects filename)

# dependencies:
	# arty.py

# input: none
# output: rosetta score of structures

##############################################################################
# calculates rosetta score of structures
##############################################################################

import os
from datetime import datetime
import argparse
import arty

dt_start = datetime.now()

parser = argparse.ArgumentParser(description = "this script calculates rosetta scores")

parser.add_argument("-t", "--target", help = "indicates that the structures are targets", action = "store_true")
parser.add_argument("-n", help = "number of structures", required = True)
parser.add_argument("-id", help = "pdb id of target", default = "1acb")
parser.add_argument("-em", help = "indicates that the structures are energy minimized (affects filename)", action = "store_true")
parser.add_argument("-r", help = "indicates that the structures are rosetta relaxed (affects filename)", action = "store_true")

args = parser.parse_args()

if args.target:

	with open("pdb_ids.txt") as input:

		for line in input:

			pdb_id = line.strip()

			if args.em:

				filename = "em_" + pdb_id + "_complex_H.pdb"

			elif args.r:

				filename = "relaxed_" + pdb_id + "_complex_H.pdb"

			else:

				filename = pdb_id + "_complex_H.pdb"

			print (filename)

			os.system("/gpfs/loomis/apps/avx/software/Rosetta/3.12/main/source/bin/score_jd2.static.linuxgccrelease -in:file:s " + filename + " -ignore_zero_occupancy false > score_output.txt")

	os.system("mv score.sc rosetta_scores_targets.txt")

else:

	for i in range(1, int(args.n) + 1):

		if args.em:

			filename = "em_" + args.id + "_model_" + str(i) + "_complex_H.pdb"

		elif args.r:

			filename = "relaxed_" + args.id + "_model_" + str(i) + "_complex_H.pdb"

		else:

			filename = args.id + "_model_" + str(i) + "_complex_H.pdb"

		print (filename)

		if os.path.exists(filename):

			os.system("/gpfs/loomis/apps/avx/software/Rosetta/3.12/main/source/bin/score_jd2.static.linuxgccrelease -in:file:s " + filename + " -ignore_zero_occupancy false > score_output.txt")

			dt = datetime.now()

	os.system("mv score.sc rosetta_scores_" + args.id + ".txt")

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y")

str_runtime = arty.get_runtime(dt_start, dt)

print ("calc_rosetta_score.py finished at " + str_dt + " (runtime " + str_runtime + ")")