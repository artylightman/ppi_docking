function [] = generate_surface(PDB_name)
PDB_name1 = strcat(PDB_name,'.txt')
position_file_name = strcat(PDB_name,'_total_position.txt')
f_write = fopen(position_file_name,'w');

f2 = fopen(strcat(PDB_name,'_param.lua'),'w');
position = strcat('positionfile = ',position_file_name);
fprintf(f2,'positionfile = ''%s''\n',position_file_name);
fprintf(f2,'readfile = ''vor_code/read.lua''\n\n');

fprintf(f_write, '#label\n');
size_number = 7;
filename = fopen(PDB_name1)
data = textscan(filename,'%s %s %d %f %f %f %f');

xyz(:,1) = data{5};
xyz(:,2) = data{6};
xyz(:,3) = data{7};
XtalPosition = xyz;
sizes_all = data{4};
sizes_square_all = sizes_all.^2;

% Load PDB
% if exist(strcat(InputFilename,PDB_name1,'_H.mat'))  
%     PDB_name = strcat( InputFilename,PDB_name1, '_H.mat');
%     PDB_name
%     load(PDB_name);
% else
%     pdbstruct = pdbread(strcat(InputFilename,PDB_name1, '_H.pdb'));
%     fprintf(strcat(InputFilename,PDB_name1, '_H.pdb\n'))
%     x = pdbstruct.Model.Atom;
%     tempModel1=struct2cell(x);
%     tempModel2=reshape(tempModel1,size(tempModel1,1),size(tempModel1,3));
%     save(strcat(InputFilename,PDB_name1,'_H.mat'), 'tempModel2');
% end
% 
% ind1 = strcmp(tempModel2(:,3),'B');
% tempModel2 = tempModel2(~ind1,:);
% tempModel2(:,1) = num2cell([1:size(tempModel2,1)]);
% 
% res_ids = cell2mat(tempModel2(:,6));
% 
% 
% %relabel residue ids of second chain so all > ids of 1st chain
% for i = 2:size(tempModel2,1)
%     if res_ids(i,1) < res_ids(i-1,1)
%         res_ids(i:size(tempModel2,1),1) = res_ids(i:size(tempModel2,1),1) + 1000;
%     end
% end
% tempModel2(:,6) = num2cell(res_ids);
% 
% tempModel2(:,1) = num2cell([1:size(tempModel2,1)]);
% 
% %Add sizes to everything
% tempModel2 = add_sizes_protein(tempModel2,size_number);
% sizes_all = cell2mat(tempModel2(:,12));
% sizes_square_all = sizes_all.^2;



% XtalPosition = cell2mat(tempModel2(:,8:10));
xmin = min(XtalPosition(:,1));
xmax = max(XtalPosition(:,1));
ymin = min(XtalPosition(:,2));
ymax = max(XtalPosition(:,2));
zmin = min(XtalPosition(:,3));
zmax = max(XtalPosition(:,3));


xcenter = (xmin+xmax)/2.0;
ycenter = (ymin+ymax)/2.0;
zcenter = (zmin+zmax)/2.0;

XtalPosition(:,1) = XtalPosition(:,1) - xcenter;
XtalPosition(:,2) = XtalPosition(:,2) - ycenter;
XtalPosition(:,3) = XtalPosition(:,3) - zcenter;

xmin = 1*min(XtalPosition(:,1))
xmax = 1*max(XtalPosition(:,1))
ymin = 1*min(XtalPosition(:,2))
ymax = 1*max(XtalPosition(:,2))
zmin = 1*min(XtalPosition(:,3))
zmax = 1*max(XtalPosition(:,3))

fprintf(f2,'xmin = %10.6f \n', xmin);
fprintf(f2,'ymin = %10.6f \n', ymin);
fprintf(f2,'zmin = %10.6f \n\n', zmin);

fprintf(f2,'xmax = %10.6f \n', xmax);
fprintf(f2,'ymax = %10.6f \n', ymax);
fprintf(f2,'zmax = %10.6f \n\n', zmax);

fprintf(f2,'epsilon = 1e-6\n\n');

fprintf(f2,'boundary = "none"\n\n');
fprintf(f2,'postprocessing = true\n');
fprintf(f2,'savesurface = false\n');
fprintf(f2,'savereduced = false\n');
fprintf(f2,'savepoly = false\n');
fclose(f2);
% tempModel2(:,8:10) = num2cell(XtalPosition);
resid = data{3};
resid_list = unique(resid);

[x,y,z] = sphere(20);

for ii = 1:length(resid_list)
    ind0 = resid == resid_list(ii);
    residue_position = XtalPosition(ind0,:);
    residue_radii = sizes_all(ind0,:);
    residue_radii_square = sizes_square_all(ind0,:);
    if size(residue_position,1) == 1
        continue;
    end
    for jj = 1:size(residue_position,1)
        
        x_atom = x*residue_radii(jj,1) + residue_position(jj,1);
        y_atom = y*residue_radii(jj,1) + residue_position(jj,2);
        z_atom = z*residue_radii(jj,1) + residue_position(jj,3);

        for xx = 1:length(x_atom)
            for yy = 1:length(x_atom)
                this_surf_dot = [x_atom(xx,yy),y_atom(xx,yy),z_atom(xx,yy)];
                
                
                keep_dot = 1;
                for ll = 1:size(residue_position,1)
                    if ll~=jj
                        dist = (this_surf_dot(1)-residue_position(ll,1))^2+ (this_surf_dot(2)-residue_position(ll,2))^2+(this_surf_dot(3)-residue_position(ll,3))^2;
                        if dist < residue_radii_square(ll)
                            keep_dot = 0;
                            break;
                        end
                    end
                end
                if keep_dot == 1
                            fprintf(f_write, ('%5d %8.3f %8.3f %8.3f\n'), resid_list(ii),this_surf_dot(1), this_surf_dot(2), this_surf_dot(3));
      
                end
            end
        end
%         for ll = 1:length(residue_position)
%             
%         end
    end
    
    
end


fclose(filename);
fclose(f_write);



end
