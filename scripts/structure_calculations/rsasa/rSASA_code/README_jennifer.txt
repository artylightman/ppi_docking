- Install naccess using : csh install.sub
- Call run_Naccess(pdb_folder, pdb_name, save_folder)
from that folder where:
pdb_folder is the full path to the folder the pdb is in (with / at the end)
pdb_name is the pdb name (ex XXXX)
save_folder is where you want to save things (also with / at end)

The code will look for the file XXXX.mat in pdb_folder.
It will output a file XXX_sasa_data.mat in save_folder
this file contains the variable each_res_data:
column 1: residue id
column 2: rSASA
column 3: SASA in protein context
column 4: dipeptide SASA

Note: residue ids in XXXX.mat must be unique across chains.
