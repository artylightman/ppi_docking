# author : artemis lightman
# date created : feb 21, 2023
# last modified: feb 21, 2023

# command line arguments
	# n - number of decoys per target
	# topn - number of top decoys to find
	# start - value to start decoy numbering from

# dependencies:
	# arty.py

# input: rosetta generated decoys
# output: moves top scored decoys into respective folders in top_decoys folder

##############################################################################
# copy and rename decoys
##############################################################################

import os
from datetime import datetime
import argparse
import arty
import glob
import pandas as pd

dt_start = datetime.now()

parser = argparse.ArgumentParser(description = "this script copies and renames the top scored rosetta generated decoys")

parser.add_argument("-n", type = int, help = "total number of decoys per target", default = 500)
parser.add_argument("-topn", type = int, help = "number of top decoys to find", default = 250)
parser.add_argument("-start", type = int, help = "value to start decoy numbering from", default = 0)

args = parser.parse_args()

a = glob.glob("*_rosetta")

pdbs = []

with open("targets.txt") as input:

	for line in input:

		arr = line.split()

		pdbs.append(arr)

for item in pdbs:

	os.mkdir("top_decoys/" + item[0])

	with open(item[0] + "_rosetta/score.sc") as input:
		
		model_nums = []
		scores = []

		for line in input:

			if line.startswith("SEQUENCE") or line.startswith("SCORE: total_score"):

				pass

			else:

				arr = line.split()
				score = float(arr[1])
				scores.append(score)
			
		for model_count in range(1, args.n + 1):
			
			if model_count < 10:
					
				model_num = "000" + str(model_count)
					
			elif model_count < 100:
				
				model_num = "00" + str(model_count)
					
			else:
					
				model_num = "0" + str(model_count)
					
			filename = item[0] + "_rosetta/" + item[0] + "_" + item[1] + "_" + model_num + ".pdb"

			if os.path.exists(filename):

				model_nums.append(model_count)

	df = pd.DataFrame({'model_number': model_nums, 'score': scores})

	df = df.sort_values(by = 'score')

	top_n = list(df.model_number[:args.topn])

	for index, val in enumerate(top_n):

		if val < 10:

			model_num = "000" + str(val)

		elif val < 100:

			model_num = "00" + str(val)

		else:

			model_num = "0" + str(val)

		filename = item[0] + "_rosetta/" + item[0] + "_" + item[1] + "_" + model_num + ".pdb"

		cmd = "cp " + filename + " top_decoys/" + item[0] + "/" + item[0] + "_model_" + str(args.start + index + 1) + ".pdb"
		os.system(cmd)

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y ")

str_runtime = arty.get_runtime(dt_start, dt)

print ("copy_decoys.py finished at " + str_dt + " (runtime " + str_runtime + ")")