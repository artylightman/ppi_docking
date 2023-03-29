# author : artemis lightman
# date created : jan 2, 2023
# last modified: jan 2, 2023

# command line arguments
    # id - pdb id of target
    # n - number of pdb structures
    # t/target - indicates that the structures are targets

# dependencies:
    # arty.py

# input: all pdb files
# output: calculates overlap energy (all data in overlap folder)

##############################################################################
# calculates overlap energy 
##############################################################################

import os
from datetime import datetime
import argparse
import arty

dt_start = datetime.now()

parser = argparse.ArgumentParser(description = "this script runs the overlap energy calculation")

parser.add_argument("-id",  help = "pdb id of target", default = "1acb")
parser.add_argument("-n", help = "number of pdb structures", type = int, default = 500)
parser.add_argument("-t", "--target", help = "indicates that the structures are targets", action = "store_true")

args = parser.parse_args()

if args.target:

	os.system("cp -r ../structure_calculations/overlap overlap")
	job_name = "overlap_targets"

else:
	
	os.system("cp -r ../../../structure_calculations/overlap overlap")
	job_name = "overlap_" + args.id

os.system("cp *_complex.txt overlap/pdb")

os.chdir("overlap/overlap_code")

os.system("python run_overlap.py")

# update submit script

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

print ("calc_overlap.py finished at " + str_dt + " (runtime " + str_runtime + ")")