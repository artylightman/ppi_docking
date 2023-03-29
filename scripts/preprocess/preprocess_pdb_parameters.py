#! /Users/jennifergaines/anaconda/bin/python
#/Library/Frameworks/Python.framework/Versions/2.7/bin/python
#import matplotlib.pyplot as plt
import Bio.PDB as pdb
import numpy as np

# Sets up the coord_size file for the protein 

import csv
with open('atypes.csv') as f:
    atypes = dict([((res,aname), atype) for res,aname,atype in csv.reader(f)])
with open('asizes_7.csv', 'r') as f:
    asizes = dict([(atype, float(sz)) for atype, sz in csv.reader(f)])
    


def get_size(atom):
	
	rname, aname = atom.parent.resname, atom.name
	if (rname, aname) not in atypes and ('X',aname) in atypes :
		rname = 'X'
	atype = atypes.get((rname, aname), None)
	return asizes.get(atype, None), rname, atype

def save_pdb_as_text(struc, size_func = get_size):
	"Save the pdb as a text file"
	model = spdb[0]
	chain = model
	#atoms = [a for a in chain.get_atoms() ]
	atoms = [a for a in chain.get_atoms() if pdb.is_aa(a.parent)]
	
	parents = []
	counter = 0;
	for  a in chain.get_atoms() :
		if pdb.is_aa(a.parent):
			parents.append(a.parent)
			counter = counter + 1;
		

	sizes = [];
	resnames = [];
	atype = [];
	for i in range (0, len(atoms)):
		a = atoms[i] 
		s,r,at = size_func(a)
		if s is None:
			s = 1.4
		sizes.append(s)
		resnames.append(r)
		atype.append(at)


	xyzs = [(a.coord) for a in atoms]
	xyzarr = np.array(xyzs)

	maxcoord = np.amax(abs(xyzarr))
	lowerlim = tuple([-maxcoord*2]*3)
	upperlim = tuple([maxcoord*2]*3)


	#f = open('/Users/jennifergaines/Documents/Protherm/all_pdb/' + file_name + '.txt', 'w')
	f = open(folder_name+file_name + '.txt', 'w')

	for i in range(0, len(sizes)) :
		if parents[i].get_id()[2] != 'A' and parents[i].get_id()[2] != 'B'  and parents[i].get_id()[2] != 'C':
			if parents[i].get_id()[2] != ' ':
				print('Bad file spot ' )
				print( parents[i].get_id())
				id_num = parents[i].get_id()[1]*10  + int(parents[i].get_id()[2])
				f.write('{:5s} {:4s} {:2d} {:3f} '.format(atoms[i].get_name(), parents[i].get_resname(), parents[i].get_id()[1] ,sizes[i]))
				f.write('{:3f} {:3f} {:3f} \n'.format(xyzarr[i][0], xyzarr[i][1], xyzarr[i][2]))
			else:
				f.write('{:5s} {:4s} {:2d} {:3f} '.format(atoms[i].get_name(), parents[i].get_resname(), parents[i].get_id()[1], sizes[i]))
				f.write('{:3f} {:3f} {:3f} \n'.format(xyzarr[i][0], xyzarr[i][1], xyzarr[i][2]))

	f.close()

############################
# Main function starts here	
############################	

import sys
print(sys.argv[1:])

hiq = open(sys.argv[2], 'r')
parser = pdb.PDBParser()
folder_name = str(sys.argv[1])

# The first structure
for hiq_files in range(0,int(sys.argv[3])):
	file_name = hiq.readline().rstrip()
	file_name = file_name.lower()
	print(folder_name+file_name+"_H.pdb")
	
	spdb = parser.get_structure("new_file", folder_name+file_name+"_H.pdb")
	model =  spdb[0]
	chain = model
	atoms = [a for a in chain.get_atoms() if pdb.is_aa(a.parent)]
	save_pdb_as_text(spdb)
hiq.close()
