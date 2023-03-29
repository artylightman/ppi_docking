%{
% function run_Naccess(pdb_folder, pdb_name, save_folder)
% 
% Input: 
% pdb_folder: full path the folder with *.mat file of the pdb. End path
name with /
% pdb_name: Name of the pdb
% save_folder: full path to folder to save results in. End path name with /
%
% Output:
% *_sasa_data.mat: contains the following columns:
%   1: Residue id
%   2: rSASA
%   3: SASA of residue in context of protein
%   4: SASA in context of dipeptide (Ca,C,O of prior amino acid and N,Ca, H
of next amino acid
% 
% Notes:
% - Naccess must be installed in this folder
% - coordinate file in pdb_name.mat must have unique residue ids for each
%   residue.
% - Set up PDB file using download_preprocess_pdb.py
%}

function run_Naccess(pdb_folder, pdb_name, save_folder)

load(strcat(pdb_folder, pdb_name, '.mat'));
res_ids = unique(cell2mat(tempModel2(:,6)));
all_ids = cell2mat(tempModel2(:,6));
each_res_data = zeros(size(res_ids,1), 4);
each_res_data(:,1) = res_ids;



%% Run full protein through naccess
dipept = tempModel2;

%Save protein to PDB file
f = fopen('this_protein.pdb', 'w');
for j = 1:size(dipept,1)
    if size(dipept{j,2},2)<4
        fprintf(f,( '%-6s%5d  %-3s%4s %s%4d%s    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n'), 'ATOM', j, dipept{j,2}, dipept{j,4}, ' ', dipept{j,6}, '', dipept{j,8}, dipept{j,9}, dipept{j,10}, dipept{j,11}, dipept{j,12}, dipept{j,14});
    else
        fprintf(f, ('%-6s%5d %-4s%4s %s%4d%s    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n'), 'ATOM', j, dipept{j,2}, dipept{j,4},' ', dipept{j,6}, '', dipept{j,8}, dipept{j,9}, dipept{j,10}, dipept{j,11}, dipept{j,12}, dipept{j,14});
    end
end
fclose(f);

%Run naccess
x = strcat('./naccess ',{' '}, 'this_protein.pdb -r vdw.radii_jennifer -p 1.4 -y  -z 0.001 > test.out');
x
system(x{:});

%Load naccess output data
f = fopen(strcat('this_protein.asa'));
data = textscan(f,'%s %f %s %s  %f %f %f %f %f %f');
all_protein = data{9};
fclose(f);

%% run dipeptide SASA
for res = 1:size(res_ids,1)
    
    %Store SASA of residue in context of protein
    ind0 = ismember(all_ids,res_ids(res));
    ind1 = ismember(each_res_data(:,1), res_ids(res));
    each_res_data(ind1,3) = sum(all_protein(ind0));
    
    %Create dipeptide
    prior = {};
    next = {};
    if res > 1
        ind1 = ismember(all_ids,res_ids(res)-1);
        ind3 = ismember(tempModel2(:,2), {'C', 'O', 'CA'});
        prior = tempModel2(ind1&ind3,:);
    end
    if res<size(res_ids,1)
        ind2 = ismember(all_ids,res_ids(res)+1);
        ind3 = ismember(tempModel2(:,2), {'CA', 'N', 'H'});
        next = tempModel2(ind2&ind3,:);
    end
    dipept = [prior;tempModel2(ind0,:);next];
    dipept(:,1) = num2cell([1:size(dipept,1)]);
    
    %Write dipeptide to file
    f = fopen('this_dipept.pdb', 'w');
    for j = 1:size(dipept,1)
        if size(dipept{j,2},2)<4
            fprintf(f,( '%-6s%5d  %-3s%4s %s%4d%s    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n'), 'ATOM', j, dipept{j,2}, dipept{j,4}, ' ', dipept{j,6}, '', dipept{j,8}, dipept{j,9}, dipept{j,10}, dipept{j,11}, dipept{j,12}, dipept{j,14});
        else
            fprintf(f, ('%-6s%5d %-4s%4s %s%4d%s    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n'), 'ATOM', j, dipept{j,2}, dipept{j,4},' ', dipept{j,6}, '', dipept{j,8}, dipept{j,9}, dipept{j,10}, dipept{j,11}, dipept{j,12}, dipept{j,14});
        end
    end
    fclose(f);
    
    %Run Naccess
    x = strcat('./naccess ',{' '}, 'this_dipept.pdb -r vdw.radii_jennifer -p 1.4 -y  -z 0.01 > test.out');
    system(x{:});
    
    %Load Naccess output and save SASA data
    f = fopen(strcat('this_dipept.asa'));
    data = textscan(f,'%s %f %s %s %f %f %f %f %f %f');
    all_new = data{9};
    ind0 = ismember(cell2mat(dipept(:,6)), res_ids(res));
    ind1 = ismember(each_res_data(:,1), res_ids(res));
    each_res_data(ind1,4) = sum(all_new(ind0));
    fclose(f);
       
    
end

%Calculate rSASA
each_res_data(:,2) = each_res_data(:,3)./each_res_data(:,4);

%Save results
save(strcat(save_folder,pdb_name, '_sasa_data.mat'), 'each_res_data');

end

