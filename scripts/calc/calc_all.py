# author : artemis lightman
# date created : jan 2, 2023
# last modified: jan 2, 2023

# command line arguments
	# id - pdb of target
	# n - number of pdb structures
	# path - path to rsasa folder
	# t/target - indicates whether structures are targets
	# em - indicates whether structures are energy minimized

# dependencies:
	# arty.py, calc_rsasa.py, calc_single_volume.py, calc_volume_voronoi.py, calc_volume_grid.py, calc_overlap.py, calc_ialign.py, packing_tasklist.sh 

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
parser.add_argument("-path", help = "path to rSASA folder", required = True)
parser.add_argument("-t", "--target", help = "indicates that the structures are targets", action = "store_true")
parser.add_argument("-em", help = "indicates that the structures are energy minimized", action = "store_true")

args = parser.parse_args()

# example path: interfaces/targets/rsasa

num_pdb = str(args.n)

if args.target:

	os.system("python calc_rsasa.py -path " + args.path + " --target")
	os.system("python calc_volume_voronoi.py -n " + num_pdb + " --target")
	os.system("python calc_overlap.py -n " + num_pdb + " --target")
	os.system("python calc_volume_grid.py -n " + num_pdb + " --target")

else:

	os.system("python calc_rsasa.py -path " + args.path + " -id " + args.id)
	os.system("python print_model_names.py")
	os.system("python calc_overlap.py -n " + num_pdb + " -id " + args.id)
	os.system("python calc_volume_grid.py -n " + num_pdb + " -id " + args.id)
	os.system("python calc_volume_voronoi.py -n " + num_pdb + " -id " + args.id)

	if args.em:

		os.system("python calc_ialign.py -em -id " + args.id)

	else:

		os.system("python calc_ialign.py -id " + args.id)

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y")

str_runtime = arty.get_runtime(dt_start, dt)

print ("calc_all.py finished at " + str_dt + " (runtime " + str_runtime + ")")