# author : artemis lightman
# date created : mar 8, 2023
# last modified: mar 9, 2023

# command line arguments
	# n - number of pdb files to preprocess
	# t/target - indicates the structures are targets

# dependencies:
	# arty.py, download_preprocess_pdb.py, bash_edge_code_script_local.sh

# input: all files with pdb extension
# output: pdb files for both chains and complex with hydrogens added, and text file version

##############################################################################
# completes pdb preprocessing
##############################################################################

import glob
import os
from datetime import datetime
import argparse
import arty
from Bio.PDB import *

dt_start = datetime.now()

parser = argparse.ArgumentParser(description = "this script completes all PDB preprocessing")

parser.add_argument("-n", type = int, help = "number of PDB files to preprocess", default = 500)
parser.add_argument("-t", "--target", help = "indicates that the structures are targets", action = "store_true")

args = parser.parse_args()

class chainSelect(Select): #filter for only items in chain specified

	def __init__(self, chain_id):

		super()

		self.chain_id = chain_id

	def accept_chain(self, chain):
        
		if chain.id == self.chain_id:
            
			return 1
        
		else:
            
			return 0

def renumber(filenames): # renumber atoms and residues to remove hydrogens and fix mislabelled residues

	for filename in filenames:

		no_ext = filename[:-4]

		output = open(no_ext + "_corrected.pdb", "w")

		res_num = 0
		res_num_prev = 0
		atom_num = 1

		with open(filename) as input:

			for line in input:

				if line.startswith("ATOM"):

					res_num_curr = line[22:27]
					atom_type = line[76:78]

					if "H" not in atom_type:

						if (res_num_curr != res_num_prev):

							res_num += 1

						if args.target:

							new_line = line[:6] + '{:>5}'.format(atom_num) + line[11:22] + '{:>4}'.format(res_num) + " " + line[27:]

						else:

							new_line = line[:6] + '{:>5}'.format(atom_num) + line[11:22] + '{:>4}'.format(res_num) + " " + line[27:54] + "  1.00  0.00" + "\n"

						atom_num += 1

						res_num_prev = res_num_curr

						output.write(new_line)

				else:

					continue

		output.close()

		os.system("rm " + filename)

		os.system("mv " + filename[:-4] + "_corrected.pdb " + filename)

def split_chains(filenames):

	for filename in filenames:

		chains = []

		structure_name  = filename[:-4]

		print ("working on pdb complex " + structure_name)

		###### pdb parsing and setup ###########

		parser = PDBParser()
		structure = parser.get_structure("structure", filename)

		for model in structure:
				
			for chain in model:
					
				chains.append(chain.id)

		# write out pdb files

		io = PDBIO()
		io.set_structure(structure)
		io.save(structure_name + "_complex.pdb")

		for i in range(len(chains)):

			io.save(structure_name + "_" + chains[i] + ".pdb", chainSelect(chain_id=chains[i]))

		os.system("rm " + filename)

def print_name(filenames):

	output = open("pdb.txt", "w")
	output_complex_only = open("pdb_complex_only.txt", "w")

	for filename in filenames: # write out filenames to be preprocessed

		if "complex" in filename:
            
			output_complex_only.write(filename[:-4] + "\n")

		output.write(filename[:-4] + "\n")

def remove_excess(filenames): # remove unnecessary pdb files

	if args.target:

		os.mkdir("noH")

	for filename in filenames:

		if args.target:

			if "ordered" in filename:

				os.system("rm " + filename)

			elif "noH" in filename and "_complex_noH" not in filename:

				os.system("mv " + filename + " noH/" + filename)

			elif "_H.pdb" not in filename:

				os.system("rm " + filename)

			elif "_H.pdb" in filename and (os.path.exists(filename[:-6] + "_H_H.pdb") == True):

				os.system("rm " + filename)


		else:

			if "noH" in filename or "ordered" in filename:

				os.system("rm " + filename)

			elif "_H.pdb" not in filename:

				os.system("rm " + filename)

			elif "_H.pdb" in filename and (os.path.exists(filename[:-6] + "_H_H.pdb") == True):

				os.system("rm " + filename)

a = glob.glob("*.pdb")

if len(a) == 0:

	print ("no pdb files found. exiting.")

else:

	renumber(a)
	split_chains(a)

	b = glob.glob("*.pdb")

	print_name(b)

	with open("pdb.txt", "r") as input: # count number of total pdb structures

		for count, line in enumerate(input):

			pass

	num_files = count + 1

	os.system("python3 download_preprocess_pdb.py ./ pdb.txt " + str(num_files) + " ./") # run reduce

	os.system("bash bash_edge_code_script_local.sh pdb_complex_only.txt ./ " + str(args.n) +  " vol_code") # generate tasklist and txt file format

	c = glob.glob("*.pdb")

	remove_excess(c)

dt = datetime.now()

str_dt = dt.strftime("%H:%M:%S on %d %B, %Y ")

str_runtime = arty.get_runtime(dt_start, dt)

print ("preprocess.py finished at " + str_dt + " (runtime " + str_runtime + ")")