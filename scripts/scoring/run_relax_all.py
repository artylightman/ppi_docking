# author : artemis lightman
# date created : feb 6, 2023
# last modified: feb 6, 2023

# command line arguments
	# name - name of pdb structure
	# id - pdb id for submitted job name

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
import glob

dt_start = datetime.now()

parser = argparse.ArgumentParser(description = "this script completes all structure calculations")

parser.add_argument("-n",  help = "number of structures to relax", required = True)

parser.add_argument("-id",  help = "pdb id for job name", required = True)

args = parser.parse_args()

output = open("tasklist_relax.sh", "w")

filenames = glob.glob("*_complex_H.pdb")

for item in filenames:

    output.write("python run_relax.py -name " + item + " \n")

job_name = "relax_" + args.id

output = open("submit_relax_corrected.sh", "w")

with open("submit_relax.sh") as input:
	
	for line in input:
		
		if line.startswith("#SBATCH --job-name="):

			line = "#SBATCH --job-name=" + job_name + "\n"
			
		elif line.startswith("#SBATCH --array=1-"):

			line = "#SBATCH --array=1-" + str(args.n) + "\n"

		
		output.write(line)
		
output.close()

os.system("rm submit_relax.sh")
os.system("mv submit_relax_corrected.sh submit_relax.sh")

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y")

str_runtime = arty.get_runtime(dt_start, dt)

print ("run_relax_all.py finished at " + str_dt + " (runtime " + str_runtime + ")")

