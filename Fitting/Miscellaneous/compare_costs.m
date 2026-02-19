clear all; close all; clc
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
[username, githubfolder] = get_paths();
iFs = [2,3,5,6,7,8,11];

mcodes = [1 1 1; 1 1 1];
v = {'','_v2'};

for j = 1:size(mcodes,1)
    mcode = mcodes(j,:);
    
for i = 1:length(iFs)
    [output_mainfolder, filename, opt_type, ~] = get_folder_and_model(mcode);
    foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iFs(i)}];

    cd(foldername)
    load(['parms_', filename, v{j}, '.mat'], 'newparms', 'optparms', 'out', 'bnds')


    J(i,j) = out.J;
end
end

%%
figure(1)
plot(J)
% legend('FD', 'Coop', 'Regular', 'Hill', 'location', 'best')