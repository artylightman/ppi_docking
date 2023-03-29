# author : artemis lightman
# date created : jan 23, 2023
# last modified: mar 29, 2023

# command line arguments
    # em - indicates whether structures are energy minimized (affects filename)
	# r - indicates whether structures are rosetta relaxed (affects filename)
	# n - number of structures
	# id - pdb id of target

# dependencies:
	# arty.py, itscsorepro

# input: none
# output: hdock score of structures

##############################################################################
# calculates hdock score of structures
##############################################################################

import os
from datetime import datetime
import argparse
import arty

dt_start = datetime.now()

parser = argparse.ArgumentParser(description = "this script calculates hdock scores")

parser.add_argument("-em", help = "indicates that the structures are energy minimized (affects filename)", action = "store_true")
parser.add_argument("-r", help = "indicates that the structures are rosetta relaxed (affects filename)", action = "store_true")
parser.add_argument("-n", help = "number of structures", required = True)
parser.add_argument("-id", help = "pdb id of target", default = "1acb")

args = parser.parse_args()

for i in range(1, int(args.n) + 1):

	if args.em:

		filename = "em_" + args.id + "_model_" + str(i) + "_complex_H.pdb"

	elif args.r:

		filename = "relaxed_" + args.id + "_model_" + str(i) + "_complex_H.pdb"

	else:

			filename = args.id + "_model_" + str(i) + "_complex_H.pdb"

	print (filename)

	if os.path.exists(filename):

		os.system("./DockQ.py " + filename + " ../../../targets/" + args.id + "_complex_H.pdb >> dockq_scores_" + args.id + ".txt")

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y")

str_runtime = arty.get_runtime(dt_start, dt)

print ("calc_dockq.py finished at " + str_dt + " (runtime " + str_runtime + ")")