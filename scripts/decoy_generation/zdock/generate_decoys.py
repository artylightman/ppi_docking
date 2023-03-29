# author : artemis lightman
# date created : jan 3, 2023
# last modified: jan 3, 2023

# command line arguments
	# n - number of decoys
	# id - pdb id of target
	# c - chains in target

# dependencies:
	# arty.py

# input: none
# output: decoys for given target

##############################################################################
# generates decoys for a given target
##############################################################################

import os
from datetime import datetime
import argparse
import arty

dt_start = datetime.now()

parser = argparse.ArgumentParser(description = "this script generates decoys for a given target using zdock")

parser.add_argument("-n", type = int, help = "number of decoys to generate", default = 500)
parser.add_argument("-id", help = "pdb id of target", required = True)
parser.add_argument("-c", help = "chains in target", required = True)

args = parser.parse_args()

chain_one = args.c[0]
chain_two = args.c[1]

num_models = args.n

receptor_pdb = "complexes/" + args.id + "_" + chain_one + "_noH.pdb"
receptor_m_pdb = args.id + "_zdock/receptor_m.pdb"
ligand_pdb = "complexes/" + args.id + "_" + chain_two + "_noH.pdb"
ligand_m_pdb = args.id + "_zdock/ligand_m.pdb"

os.mkdir(args.id + "_zdock")

os.system("./mark_sur " + receptor_pdb + " " + receptor_m_pdb)
os.system("./mark_sur " + ligand_pdb + " " + ligand_m_pdb)
os.system("cp create.pl " + args.id + "_zdock/create.pl")
os.system("cp create_lig " + args.id + "_zdock/create_lig")
os.system("cp zdock " + args.id + "_zdock/zdock")

os.chdir(args.id + "_zdock")

os.system("./zdock -R receptor_m.pdb -L ligand_m.pdb -o zdock.out -D -N " + str(num_models))

os.system("perl create.pl zdock.out")

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y ")

str_runtime = arty.get_runtime(dt_start, dt)

print ("generate_decoys.py finished at " + str_dt + " (runtime " + str_runtime + ")")