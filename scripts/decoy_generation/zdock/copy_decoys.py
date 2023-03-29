# author : artemis lightman
# date created : jan 3, 2023
# last modified: jan 3, 2023

# command line arguments
	# n - number of decoys per target

# dependencies:
	# arty.py

# input: zdock generated decoys
# output: copies decoys to respective folders in all_decoys directory

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

a = glob.glob("*_zdock")

for item in a:
	os.mkdir("all_decoys/" + item[0:4])
	for i in range(1, (args.n + 1)):
		filename = item + "/complex." + str(i) + ".pdb"

		if os.path.exists(filename):
			cmd = "cp " + filename + " all_decoys/" + item[0:4] + "/" + item[0:4] + "_model_" + str(i) + ".pdb"
			os.system(cmd)

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y ")

str_runtime = arty.get_runtime(dt_start, dt)

print ("copy_decoys.py finished at " + str_dt + " (runtime " + str_runtime + ")")