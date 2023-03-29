# author : artemis lightman
# date created : jan 2, 2023
# last modified: jan 2, 2023

# command line arguments
    # id - pdb id of target
    # n - number of pdb structures to process
    # t/target - indicates whether structures are targets

# dependencies:
    # arty.py

# input: pdb files in text form
# output: calculates residue volume using grid method (in packing folder)

##############################################################################
# calculates residue volumes in structures using grid method
##############################################################################

import sys
import glob
import os
from datetime import datetime
import argparse
import arty

dt_start = datetime.now()

parser = argparse.ArgumentParser(description = "this script calculates volume using the grid method")

parser.add_argument("-id",  help = "pdb id of target", default = "1acb")

parser.add_argument("-n", help = "number of pdb structures", required = True)

parser.add_argument("-t", "--target", help = "indicates that the structures are targets", action = "store_true")

args = parser.parse_args()

# copy volume files to desired folder

if args.target:

	os.system("cp -r ../structure_calculations/packing packing")
	job_name = "grid_targets"

else:

	os.system("cp -r ../../../structure_calculations/packing packing")
	job_name = "grid_" + args.id

os.system("cp *_complex.txt packing/vol_code")

os.system("cp packing_tasklist.sh packing/tasklist.sh")

os.chdir("packing/vol_code")

os.system("chmod 777 *")

os.chdir("..")

# update submit script to have correct job name and number of structures

output = open("submit_corrected.sh", "w")

with open("submit.sh") as input:
	
	for line in input:
		
		if line.startswith("#SBATCH --job-name="):

			line = "#SBATCH --job-name=" + job_name + "\n"
			
		elif line.startswith("#SBATCH --array=1-"):
			
			line = "#SBATCH --array=1-" + str(args.n) + "\n"

		
		output.write(line)
		
output.close()

os.system("rm submit.sh")
os.system("mv submit_corrected.sh submit.sh")

os.system("sbatch submit.sh")

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y")

str_runtime = arty.get_runtime(dt_start, dt)

print ("calc_volume_grid.py finished at " + str_dt + " (runtime " + str_runtime + ")")