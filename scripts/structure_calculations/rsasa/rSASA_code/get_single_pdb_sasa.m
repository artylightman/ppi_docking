function get_single_pdb_sasa(mdir,pdbdir,pdbfname,savedir)

% protein directories (multiple chains)
pdbf = [pdbdir pdbfname];
pdb_name = pdbfname(1:end-4);

% get pdb struct
fprintf('Loading pdb...');
pdbstruct = pdbread(pdbf);
fprintf('pdb loaded.\n');

% turn into cell
data = pdbstruct.Model.Atom;
tempModel1 = struct2cell(data);
tempModel2 = reshape(tempModel1, size(tempModel1,1), size(tempModel1,3))';
save([mdir pdb_name '.mat'],'tempModel2');

% run naccess
fprintf('Running function run_Naccess on protein...');
run_Naccess(mdir, pdb_name, savedir);
fprintf('function completed.\n');

end
