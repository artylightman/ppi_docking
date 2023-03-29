# Bash script for running edge code 
#bash bash_edge_code_script_local.sh pdb_list.txt folder1 X c1_folder c2_folder

#pdb_list.txt: list of PDB codes, in same folder as this code
#folder1: folder containing PDBs
#X: replace by integer stating number of PDBS in PDB file list
#c1_folder: full path to cluster folder that will contain volume code
#c2_folder: full path to cluster folder that will contain PDB data
# Input :
# 1. name of a file with a list of pdb names. must be full path name
# 2. folder to save everything to
# 3. Number of pdbs to run

# Ouput:
# .txt file for each protein
# Task list for ./vor tasks

# Note:
# All pdb files must already be in the folder
# You must run download_preprocess_pdb.py first to remove heteroatoms and add hydrogens to the proteins

#!/bin/bash 
rm packing_tasklist.sh
rm process.sh

file1=$1 		#File containing list of PBDB codes
file2="packing_tasklist.sh"
folder_name=$2		#Folder containing PDB files
num_pdb=$3		#Number of PDB to run
cluster_folder=$4	#Folder on cluster where code will be stored
cluster_data=$5		#Folder on cluster where PDBs and output will be
file3="process.sh"

file_ending4=".txt"
space=" "

# Create text files for the proteins
python3 preprocess_pdb_parameters.py $folder_name $file1 $num_pdb

# Make task list for getting atomic volumes and surface/core designation
while IFS='' read -r line || [[ -n "$line" ]];
	do
		run_loop=1
		echo $line
		line=$(echo "$line" | tr '[:upper:]' '[:lower:]')

		a=($(wc $folder_name$line$file_ending4))
		num_lines=${a[0]}
		
		while [ $run_loop -le 1 ]
			do
				echo -n "source ~/.bashrc; cd " >> $file2
				echo -n $cluster_folder >> $file2
				echo -n "; ./vor " >> $file2
				echo -n $cluster_data >> $file2
				echo -n $line >> $file2
				echo -n "$space" >> $file2
				echo -n $num_lines >> $file2
				echo -n "$space" >> $file2
				echo -n $run_loop >> $file2
				echo ";" >> $file2
				run_loop=$((run_loop+1))
			done
	done<$file1

#Create processing tasks
while IFS='' read -r line || [[ -n "$line" ]];
	do
		run_loop=1
		echo $line
		line=$(echo "$line" | tr '[:upper:]' '[:lower:]')

		a=($(wc $folder_name$line$file_ending4))
		num_lines=${a[0]}
		
		echo -n "source ~/.bashrc; cd " >> $file3
		echo -n $cluster_folder >> $file3
		echo -n '; module load MATLAB/2016b; matlab -nosplash -nodisplay -nojvm -r "process_volume_output(' >> $file3
		echo -n "'" >> $file3
		echo -n $line >> $file3
		echo -n "'" >> $file3
		echo -n ",'" >> $file3
		echo -n $cluster_data >> $file3
		echo -n "'," >> $file3
		echo -n $num_lines >> $file3
		echo -n ')"' >> $file3
		echo ";" >> $file3
		run_loop=$((run_loop+1))

	done<$file1

