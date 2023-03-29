import glob

files = glob.glob('../pdb/*.txt')

with open('tasklist.sh', 'w') as f:
	for file in files:
		f.write('python calc_single_overlap_energy.py '+file[7:]+';')
		f.write('\n')
