clear all; close all; clc
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
[username, githubfolder] = get_paths();
iFs = [2,3,5,6,7,8,11];

mcode = [1 1 1];
nv = nan(length(iFs), 10);

for i = 1:length(iFs)
    [output_mainfolder, filename, opt_type, ~] = get_folder_and_model(mcode);
    foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iFs(i)}];

    cd(foldername)
    load(['parms_', filename, '_v2.mat'], 'newparms', 'optparms', 'out', 'bnds')

    for j = 1:length(optparms)
        nv(i,j) = newparms.(optparms{j});
    end
    
    figure(1)
    nexttile
    bar(categorical(optparms), out.s)

end

%%
if ishandle(2), close(2); end
figure(2)
for j = 1:length(optparms)
    nexttile
    bar(iFs, nv(:,j));
    title(optparms{j})
    
    yline(bnds.(optparms{j})(1), 'k--')
    yline(bnds.(optparms{j})(2), 'k--')
end

%% quantify variability
std(nv)./mean(nv)