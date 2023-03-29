%% Get sasa for dataset of pdb files

clear;
close all;
clc;


% get pdb dir
parentdir = ['/home/gm633/project/capri_new/decoys/decoys_capri_30/rSASA/'];
pdbdir = [parentdir 'pdb/'];
% pdbdir = [parentdir 'nmr_separate_structures/'];

% get savedir
savedir = [parentdir 'rSASA_data/'];
if ~exist(savedir,'dir')
    mkdir(savedir);
end
% get mdir
mdir = [parentdir 'pdb_mat/'];

% get sasa for all proteins in pdb
flist = dir(fullfile([pdbdir '*_H.pdb']));
NF = length(flist);

for ff = 1:NF
    fprintf('on pdb %s\n',flist(ff).name);
    get_single_pdb_sasa(mdir,pdbdir,flist(ff).name,savedir);
end
