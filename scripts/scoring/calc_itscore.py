# author : artemis lightman
# date created : jan 2, 2023
# last modified: mar 9, 2023

# command line arguments
	# t/target - indicates whether structures are targets
	# em - indicates whether structures are energy minimized (affects filename)
	# r - indicates whether structures are rosetta relaxed (affects filename)
	# n - number of structures
	# id - pdb id of target

# dependencies:
	# arty.py, itscsorepro

# input: none
# output: hdock/itscorepro score of structures

##############################################################################
# calculates hdock/itscorepro score of structures
##############################################################################

import os
from datetime import datetime
import argparse
import arty

dt_start = datetime.now()

parser = argparse.ArgumentParser(description = "this script calculates hdock scores")

parser.add_argument("-t", "--target", help = "indicates that the structures are targets", action = "store_true")
parser.add_argument("-em", help = "indicates that the structures are energy minimized (affects filename)", action = "store_true")
parser.add_argument("-r", help = "indicates that the structures are rosetta relaxed (affects filename)", action = "store_true")
parser.add_argument("-n", help = "number of structures", required = True)
parser.add_argument("-id", help = "pdb id of target", default = "1acb")

args = parser.parse_args()

def reformat_coords(filename):

	output = open(filename + "_mod.pdb", "w")

	with open(filename + ".pdb") as input:
		
		for line in input:
			
			if line.startswith("ATOM"):
			
				pdb_f_start = line[:30]
				
				x_coord = line[30:38]
				
				y_coord = line[38:46]
				
				z_coord = line[46:54]
							
				new_data = pdb_f_start + x_coord + " " + y_coord + " " + z_coord
				
				output.write(new_data + "\n")

	output.close()

if args.target:

	with open("pdb_ids.txt") as input:

		for line in input:

			pdb_id = line.strip()

			if args.em:

				reformat_coords("em_" + pdb_id + "_complex_H")

				filename = "em_" + pdb_id + "_complex_H_mod.pdb"

			elif args.r:
				
				reformat_coords("relaxed_" + pdb_id + "_complex_H")

				filename = "relaxed_" + pdb_id + "_complex_H_mod.pdb"

			else:

				reformat_coords(pdb_id + "_complex_H")

				filename = pdb_id + "_complex_H_mod.pdb"

			print (filename)

			os.system("./ITScorePro " + filename + " >> hdock_scores_targets.txt")

else:

	for i in range(1, int(args.n) + 1):

		if args.em:

			reformat_coords("em_" + args.id + "_model_" + str(i) + "_complex_H")

			filename = "em_" + args.id + "_model_" + str(i) + "_complex_H_mod.pdb"

		elif args.r:

			reformat_coords("relaxed_" + args.id + "_model_" + str(i) + "_complex_H")

			filename = "relaxed_" + args.id + "_model_" + str(i) + "_complex_H_mod.pdb"

		else:

			reformat_coords(args.id + "_model_" + str(i) + "_complex_H")

			filename = args.id + "_model_" + str(i) + "_complex_H_mod.pdb"

		print (filename)

		if os.path.exists(filename):

			os.system("./ITScorePro " + filename + " >> hdock_scores_" + args.id + ".txt")

os.system("rm -rf *_mod.pdb")

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y")

str_runtime = arty.get_runtime(dt_start, dt)

print ("calc_itscore.py finished at " + str_dt + " (runtime " + str_runtime + ")")