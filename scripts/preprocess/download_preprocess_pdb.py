##############################################################
# download_preprocess_pdb.py $folder_name $file1 $num_pdb $reduce_folder
# Input: 
#  folder_name: folder with pdb files, end with /
#  file1: full path to file containing list of PDB names
#  num_pdb: number of PDBs to run
#  reduce_folder: full path to folder containing reduce
#
# Output Files:
#  XXXX_ordered.pdb for each PDB file: PDB file with all heteroatoms, disordered atoms removed
#  XXXX_noH.pdb: PDB file with no H
#  XXXX_H.pdb: PDB file with hydrogen atoms added using reduce
# 
# Notes:
# - All PDB files must be in folder_name
##############################################################

#! /Users/jennifergaines/anaconda/bin/python
from Bio.PDB import*
import Bio.PDB as pdb
import numpy as np
from subprocess import call
import os.path
import sys


class NotDisordered(Select):
	def accept_atom(self, atom):
		if (not atom.is_disordered()) or atom.get_altloc() == 'A':
			atom.set_altloc(' ')  # Eliminate alt location ID before output.

			return True
		else:  # Alt location was not one to be output.
			return False


class NonHetSelect(Select):
    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0		


hiq = open(sys.argv[2], 'r')
		
parser=PDBParser()
folder_name = str(sys.argv[1])
reduce_folder = str(sys.argv[4])
# The first structure
for hiq_files in range(0,int(sys.argv[3])):


	file_name = hiq.readline().rstrip()		
	s = parser.get_structure("my_pdb", folder_name + file_name + ".pdb")
	io=PDBIO()
	io.set_structure(s)
	
	# Remove disordered atoms, 
	io.save(folder_name + file_name + "_ordered2.pdb", select=NotDisordered())
	s = parser.get_structure("my_pdb", folder_name + file_name + "_ordered2.pdb")
	io = PDBIO()
	io.set_structure(s)
	
	#Remove heteroatoms
	io.save(folder_name + file_name + "_ordered1.pdb",  select=NonHetSelect())
	s = parser.get_structure("my_pdb", folder_name + file_name + "_ordered1.pdb")
	model = s[0]
	chain = model
	atoms = [a for a in chain.get_atoms() if pdb.is_aa(a.parent)]

	# Renumber residues to be sequential 
	parents = []
	counter = 0;
	for  a in chain.get_atoms() :
		if pdb.is_aa(a.parent):
			parents.append(a.parent)
			counter = counter + 1;
	xyzs = [(a.coord) for a in atoms]
	xyzarr = np.array(xyzs)
	f = open(folder_name+file_name + '_ordered.pdb', 'w')
	

	# Write to PDB file
	print('ordered')
	id_counter = 0
	old_res = -1
	for i in range (0, len(atoms)):
		new_res = parents[i].get_id()[1];
		if new_res != old_res:
			id_counter = id_counter+1
		count = 0
		for atom in parents[i].get_atoms():
			count += 1
		if count > 3:
			if len(atoms[i].get_name())<4:
				f.write('{:6s}{:5d}  {:<4}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s} \n'.format('ATOM', i, atoms[i].get_name(), parents[i].get_resname(),atoms[i].get_full_id()[2],  id_counter, '',xyzarr[i][0], xyzarr[i][1], xyzarr[i][2], atoms[i].get_occupancy(), atoms[i].get_bfactor(),atoms[i].get_name()[0] ))
			else:
				f.write('{:6s}{:5d} {:<4} {:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s} \n'.format('ATOM', i, atoms[i].get_name(), parents[i].get_resname(),atoms[i].get_full_id()[2],  id_counter, '',xyzarr[i][0], xyzarr[i][1], xyzarr[i][2], atoms[i].get_occupancy(), atoms[i].get_bfactor(),atoms[i].get_name()[0] ))	
		old_res = new_res			
	f.close()
	
	# Call reduce to add hydrogen atoms
	reduce1 = reduce_folder + "./reduce -Trim -quiet " +  folder_name + file_name + "_ordered.pdb>"+ folder_name + file_name + "_noH.pdb"
	reduce2 = reduce_folder + "./reduce -quiet " +  folder_name + file_name + "_noH.pdb>" +folder_name + file_name + "_H.pdb"
	os.system(reduce1)
	os.system(reduce2)
	
hiq.close()
