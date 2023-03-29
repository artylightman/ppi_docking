# author : artemis lightman
# date created : jan 2, 2023
# last modified: feb 14, 2023

# command line arguments
	# id - pdb of target
	# n - number of pdb structures
	# em - indicates that structures are energy minimized
	# r - indicates that structures are rosetta relaxed

# dependencies:
	# arty.py 

# input: pdb structures and packing tasklist
# output: runs calculations for packing fraction, rsasa, overlap energy, and ialign (decoys only)

##############################################################################
# completes structure calculations
##############################################################################

import os
from datetime import datetime
import argparse
import arty

dt_start = datetime.now()

parser = argparse.ArgumentParser(description = "this script completes all structure calculations")

parser.add_argument("-id",  help = "pdb id of target", default = "1acb")
parser.add_argument("-n", help = "number of PDB structures", type = int, default = 500)
parser.add_argument("-em", help = "indicates that the structures are energy minimized (affects filename)", action = "store_true")
parser.add_argument("-r", help = "indicates that the structures are rosetta relaxed (affects filename)", action = "store_true")

args = parser.parse_args()

os.system("cp ../../structure_calculations/scoring/calc_itscore.py " + args.id + "/calc_itscore.py")

os.system("cp ../../structure_calculations/scoring/calc_zrank.py " + args.id + "/calc_zrank.py")

os.system("cp ../../structure_calculations/scoring/calc_rosetta_score.py " + args.id + "/calc_rosetta_score.py")

os.system("cp ../../structure_calculations/scoring/ITScorePro " + args.id + "/ITScorePro")

os.system("cp ../../structure_calculations/scoring/zrank " + args.id + "/zrank")

os.system("cp ../../structure_calculations/scoring/hdock_reformat_coords.py " + args.id + "/hdock_reformat_coords.py")

os.chdir(args.id)

if args.em:

	os.system("python calc_itscore.py -n " + str(args.n) + " -id " + args.id + " -em")

	os.system("python calc_rosetta_score.py -n " + str(args.n) + " -id " + args.id + " -em")

	os.system("python calc_zrank.py -n " + str(args.n) + " -id " + args.id + " -em")

elif args.r:

	os.system("python calc_itscore.py -n " + str(args.n) + " -id " + args.id + " -r")

	os.system("python calc_rosetta_score.py -n " + str(args.n) + " -id " + args.id + " -r")

	os.system("python calc_zrank.py -n " + str(args.n) + " -id " + args.id + " -r")

else:

	os.system("python calc_itscore.py -n " + str(args.n) + " -id " + args.id)

	os.system("python calc_rosetta_score.py -n " + str(args.n) + " -id " + args.id)

	os.system("python calc_zrank.py -n " + str(args.n) + " -id " + args.id)

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y")

str_runtime = arty.get_runtime(dt_start, dt)

print ("run_score_target.py finished at " + str_dt + " (runtime " + str_runtime + ")")

