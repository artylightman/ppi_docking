#!/usr/bin/env bash

# naccess pdb_file -p probe_size (1.4 default) -r vdw_radii_file -s std_data_file -z zslice (0.05) -[hwyfaclqb]"
# -y: keep hydrogens
# -c: contact surface
# -a: excludes summary.rsa file

# -h: keep HETATMs
# -w: keep waters
# -f: keep occupancies and B-factors in output file
# -l: long .rsa format
# -q: prints the options list
# -b: dont consider alpha carbon atoms as part of the side chain ()

./naccess test.pdb -p 1.4 -z 0.01 -y -a 


# These are matlab commands I run:
# cd(datasetDirectory);
# system('mkdir asa log');
# system(['cp {' naccessDirectory 'naccess,' naccessDirectory 'vdw.radii} ' datasetDirectory]);
# probeSize = '1.4'
# 
# cd([datasetDirectory 'PDBH_Files']);
# files = io.getFiles('pdbModH');
# cd ..;
# 
# for iFile = 1:length(files)
#     iFile
#     system(['./naccess PDBH_Files/' files{iFile} ' -p ' probeSize ' -z 0.001 -y -a']);
#     if mod(iFile,100)==0
#         system('mv *.log log/; mv *.pdb pdb/; mv *.asa asa/');
#     end
# end
# system('mv *.log log/; mv *.pdb pdb/; mv *.asa asa/; rm naccess vdw.radii');
