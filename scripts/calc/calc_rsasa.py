# author : artemis lightman
# date created : jan 2, 2023
# last modified: mar 15, 2023

# command line arguments
	# id - pdb id of target
	# t/target - indicates whether structures are targets
	# path - path to rsasa folder

# dependencies:
	# arty.py

# input: all files with pdb extension
# output: rsasa calculations (in rsasa folder)

##############################################################################
# calculates rsasa
##############################################################################

import os
from datetime import datetime
import argparse
import arty

dt_start = datetime.now()

parser = argparse.ArgumentParser(description = "this script calculates rSASA")

parser.add_argument("-id",  help = "pdb id of target", default = "1acb")
parser.add_argument("-t", "--target", help = "indicates that the structures are targets", action = "store_true")
parser.add_argument("-path", help = "path to rSASA folder (starting from project directory)", required = True)

args = parser.parse_args()

if args.target:

	os.system("cp -r ../structure_calculations/rsasa rsasa") # copy rsasa files to desired folder

else:

	os.system("cp -r ../../../structure_calculations/rsasa rsasa") # path is slightly different for models

os.system("cp *.pdb rsasa/pdb") # copy all pdb files to pdb folder within rsasa

if args.target:

	path = args.path
	job_name = "rsasa_targets"

else:

	path = args.path
	job_name = "rsasa_" + args.id

#### update database_sasa.m with correct path #########

output = open("rsasa/rSASA_code/database_sasa_corrected.m", "w")

with open("rsasa/rSASA_code/database_sasa.m") as input:

	for line in input:

		if line.startswith("parentdir = "):

			line = "parentdir = ['/gpfs/gibbs/pi/ohern/gm633/" + path + "/'];\n"

		output.write(line)

os.system("rm rsasa/rSASA_code/database_sasa.m")
os.system("mv rsasa/rSASA_code/database_sasa_corrected.m rsasa/rSASA_code/database_sasa.m")

output.close()

# update job name in submit script

submit_output = open("rsasa/submit_corrected.sh", "w")

with open("rsasa/submit.sh") as submit_input:
	
	for line in submit_input:

		if line.startswith("#SBATCH --job-name="):
			
			line = "#SBATCH --job-name=" + job_name + "\n"
			
		submit_output.write(line)
		
submit_output.close()

os.system("rm rsasa/submit.sh")
os.system("mv rsasa/submit_corrected.sh rsasa/submit.sh")

os.chdir("rsasa/rSASA_code")

os.system("csh install.sub")

os.chdir("../")

os.system("sbatch submit.sh")

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y")

str_runtime = arty.get_runtime(dt_start, dt)

print ("calc_rsasa.py finished at " + str_dt + " (runtime " + str_runtime + ")")