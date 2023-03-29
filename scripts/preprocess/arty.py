# author : artemis lightman
# date created : jan 2, 2023
# last modified: jan 2, 2023

# command line arguments
	# none

# dependencies:
	# none

# input: none
# output: none

##############################################################################
# assorted useful functions
##############################################################################

from datetime import datetime
import math
import os
import numpy as np

def get_runtime(start, end): # converts start/end datetime object to runtime in h/m/s

	runtime = end - start

	rt_seconds = runtime.total_seconds()

	rt_hrs = math.floor(rt_seconds / 3600)

	rt_min = math.floor((rt_seconds % 3600) / 60)

	rt_sec = rt_seconds % 60

	str_runtime = str(rt_hrs) + "h " + str(rt_min) + "m " + str(rt_sec) + "s"

	return str_runtime

def plot_densities(values, start, stop, num_bins): # creates probability distribution plot

	plot_bins = np.linspace(start, stop, num_bins)

	density = []

	for i in range(0, len(plot_bins)-1):

		lower = plot_bins[i]
		upper = plot_bins[i+1]

		count = 0

		for j in range(0, len(values)):

			data_element = values[j]

			if lower <= data_element < upper:

				count += 1

		density.append(count/len(values) * 100)

	return (plot_bins, density)

class Core:

	def __init__(self, residues):

		self.residues = residues

class Residue:

	def __init__(self, res_id, res_type, packing_fraction, rsasa, del_rsasa, overlap):
        
		self.res_id = res_id
		self.res_type = res_type
		self.packing_fraction = packing_fraction
		self.rsasa = rsasa
		self.del_rsasa = del_rsasa
		self.overlap = overlap

def read_feature_data(filename): # returns residue feature information for given pdb structure

	residues = []

	with open(filename) as input:

		for line in input:

			if line.startswith("Res"):

				pass

			else:

				arr = line.split()
				res_id = int(arr[0])
				res_type = arr[1]
				packing_fraction = float(arr[2])
				rsasa = float(arr[3])
				del_rsasa = float(arr[4])
				overlap = float(arr[5])
				res = Residue(res_id, res_type, packing_fraction, rsasa, del_rsasa, overlap)
				residues.append(res)

	return residues

hydrophobicity_dict = {'ARG':0, 'ASP':0.09, 'GLU':0.16, 'LYS':0.16, 'ASN':0.25,
						'GLN':0.29, 'PRO':0.39, 'HIS':0.4, 'SER':0.42, 'THR':0.48,
						'GLY':0.52, 'TYR':0.64, 'ALA':0.67, 'CYS':0.74, 'MET':0.84,
						'TRP':0.85, 'VAL':0.89, 'PHE':0.96, 'LEU':0.97, 'ILE':1.0, 'GLX': 0.29}

def get_feature_data(filenames): # calculates average feature values from raw data

	fc_list = []
	avg_packing_list = []
	avg_hydro_list = []
	overlap_list = []
	all_packing = []

	for struct in filenames:

		if os.path.exists(struct):

			residues = read_feature_data(struct)

			core_residues = [res for res in residues if res.rsasa <= 1E-2]

			core_interface_residues = [res for res in core_residues if res.del_rsasa >= 1E-3]

			interface_residues = [res for res in residues if res.del_rsasa >= 1E-3]

			if len(interface_residues) > 0:

				fraction_core = len(core_interface_residues)/len(interface_residues)

				fc_list.append(fraction_core)

			else:

				fc_list.append(0)

			total_packing = 0
			total_hydrophobicity = 0

			for res in core_interface_residues:

				overlap_list.append(res.overlap)
				total_packing += res.packing_fraction
				total_hydrophobicity += hydrophobicity_dict[res.res_type]
				all_packing.append(res.packing_fraction)

			if len(core_interface_residues) > 0:

				average_packing = total_packing / len(core_interface_residues)
				average_hydrophobicity = total_hydrophobicity / len(core_interface_residues)

				avg_packing_list.append(average_packing)
				avg_hydro_list.append(average_hydrophobicity)

			else:

				avg_packing_list.append(0)
				avg_hydro_list.append(0)

	return fc_list, avg_packing_list, avg_hydro_list, overlap_list, all_packing