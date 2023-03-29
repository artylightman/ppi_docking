# author : artemis lightman
# date created : jan 2, 2022
# last modified: jan 2, 2022

# command line arguments
	# id - pdb id of target
	# n - number of pdb structures to process
	# t/target - indicates whether structures are targets

# dependencies:
	# arty.py, calc_single_volume.py

# input: pdb files in txt format
# output: calculates voronoi volume of residues

##############################################################################
# completes voronoi volume calculation
##############################################################################

import os
from datetime import datetime
import argparse
import arty

dt_start = datetime.now()

parser = argparse.ArgumentParser(description = "this script calculates the voronoi volume of residues in structures")

parser.add_argument("-id",  help = "pdb id of target", default = "1acb")
parser.add_argument("-n", help = "number of PDB structures", type = int, default = 500)
parser.add_argument("-t", "--target", help  = "indicates that the structures are targets", action = "store_true")

args = parser.parse_args()

# copy volume/pomelo files

if args.target:

	os.system("cp -r ../structure_calculations/voronoi pomelo_vor")
	os.system("cp -r ../structure_calculations/voronoi pomelo_vor_10")
	job_name = "volume_targets"
	filename = "pdb_ids.txt"

else:

	os.system("cp -r ../../../structure_calculations/voronoi pomelo_vor")
	os.system("cp -r ../../../structure_calculations/voronoi pomelo_vor_10")
	job_name = "volume_" + args.id
	filename = "models.txt"

# copy pdb files

os.system("cp *_complex.txt pomelo_vor/generating_surfaces")
os.system("cp *_complex.txt pomelo_vor_10/generating_surfaces")

##### write tasklist #########

output = open("tasklist.sh", "w")

with open(filename) as input:
	
	for line in input:
		
		pdb_name = line.strip()
		
		output.write("python calc_single_volume.py -fname " + pdb_name + ";\n")
		
output.close()

##### update submit script to have the correct job name and number of structures #####

output = open("submit_volume_corrected.sh", "w")

with open("submit_volume.sh") as input:
	
	for line in input:
		
		if line.startswith("#SBATCH --job-name="):

			line = "#SBATCH --job-name=" + job_name + "\n"
			
		elif line.startswith("#SBATCH --array=1-"):

			line = "#SBATCH --array=1-" + str(args.n) + "\n"

		
		output.write(line)
		
output.close()

os.system("rm submit_volume.sh")
os.system("mv submit_volume_corrected.sh submit_volume.sh")

####### increase box size to identify closed voronoi cells #######

os.chdir("pomelo_vor_10")

output = open("generate_surface_10.m", "w")

with open("generate_surface.m") as input:

	for line in input:
		
		if line.startswith("xmin = 1*min(XtalPosition(:,1))"):

			line = "xmin = 1*min(XtalPosition(:,1)) - 10\n"

		elif line.startswith("xmax = 1*max(XtalPosition(:,1))"):

			line = "xmax = 1*max(XtalPosition(:,1)) + 10\n"

		elif line.startswith("ymin = 1*min(XtalPosition(:,2))"):

			line = "ymin = 1*min(XtalPosition(:,2)) - 10\n"

		elif line.startswith("ymax = 1*max(XtalPosition(:,2))"):

			line = "ymax = 1*max(XtalPosition(:,2)) + 10\n"

		elif line.startswith("zmin = 1*min(XtalPosition(:,3))"):

			line = "zmin = 1*min(XtalPosition(:,3)) - 10\n"

		elif line.startswith("zmax = 1*max(XtalPosition(:,3))"):
			
			line = "zmax = 1*max(XtalPosition(:,3)) + 10\n"
			
		output.write(line)
		
output.close()

os.system("rm generate_surface.m")

os.system("mv generate_surface_10.m generate_surface.m")

os.chdir("../")

os.system("sbatch submit_volume.sh")

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y")

str_runtime = arty.get_runtime(dt_start, dt)

print ("calc_volume_voronoi.py finished at " + str_dt + " (runtime " + str_runtime + ")")