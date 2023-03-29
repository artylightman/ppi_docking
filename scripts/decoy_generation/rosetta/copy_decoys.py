# author : artemis lightman
# date created : jan 3, 2023
# last modified: jan 3, 2023

# command line arguments
    # n - number of decoys per target

# dependencies:
    # arty.py

# input: rosetta generated decoys
# output: moves decoys into respective folders in all_decoys folder

##############################################################################
# copy and rename decoys
##############################################################################

import os
from datetime import datetime
import argparse
import arty
import glob

dt_start = datetime.now()

parser = argparse.ArgumentParser(description = "this script copies and renames generated decoys")

parser.add_argument("-n", type = int, help = "number of decoys per target", default = 500)

args = parser.parse_args()

a = glob.glob("*_rosetta")

pdbs = []

with open("targets.txt") as input:

	for line in input:

		arr = line.split()

		pdbs.append(arr)

for item in pdbs:

	os.mkdir("all_decoys/" + item[0])

	for i in range(1, (args.n + 1)):

		if i < 10:
			
			model_num = "000" + str(i)
			
		elif i < 100:

			model_num = "00" + str(i)

		else:

			model_num = "0" + str(i)

		filename = item[0] + "_rosetta/" + item[0] + "_" + item[1] + "_" + model_num + ".pdb"

		if os.path.exists(filename):

			cmd = "cp " + filename + " all_decoys/" + item[0] + "/" + item[0] + "_model_" + str(i) + ".pdb"
			os.system(cmd)

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y ")

str_runtime = arty.get_runtime(dt_start, dt)

print ("copy_decoys.py finished at " + str_dt + " (runtime " + str_runtime + ")")