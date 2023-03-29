# author : artemis lightman
# date created : feb 6, 2023
# last modified: feb 6, 2023

# command line arguments
	# name - name of pdb structure

# dependencies:
	# arty.py

# input: pdb structures to relax
# output: relaxed pdb structures

##############################################################################
# runs rosetta relax on pdb structures
##############################################################################

import os
from datetime import datetime
import argparse
import arty

dt_start = datetime.now()

parser = argparse.ArgumentParser(description = "this script completes all structure calculations")

parser.add_argument("-name",  help = "name of pdb structure", required = True)

args = parser.parse_args()

os.system("/gpfs/loomis/apps/avx/software/Rosetta/3.12/main/source/bin/relax.static.linuxgccrelease -s " + args.name + " -relax:constrain_relax_to_start_coords > relax_output.txt")

pdb_name = args.name[:-4]

os.system("mv " + pdb_name + "_0001.pdb relaxed_" + pdb_name[:-10] + ".pdb")

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y")

str_runtime = arty.get_runtime(dt_start, dt)

print ("run_relax.py finished at " + str_dt + " (runtime " + str_runtime + ")")

