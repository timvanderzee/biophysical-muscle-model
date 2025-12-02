clear all; close all; clc

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
[username, githubfolder] = get_paths();
iFs = [2,3,5,6,7,8,11];

mcode = [1 1 1];
nv = nan(length(iFs), 10);

colors = parula(length(iFs));

for i = 1:length(iFs)
    [output_mainfolder, filename, opt_type, ~] = get_folder_and_model(mcode);
    foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iFs(i)}];

    cd(foldername)
    load(['parms_', filename, '_v3.mat'], 'newparms', 'optparms', 'out', 'bnds')

    newparms.kse0 = max(newparms.kse0, .1);
    
    Fs = linspace(0,2,100);
    Ls = newparms.Lse_func(Fs, newparms);

    plot(Ls*newparms.h*1e9, Fs, 'color', colors(i,:), 'DisplayName', num2str(iFs(i))); hold on
    
    kse(i) = newparms.kse;
    kse0(i) = newparms.kse0;
    
end

yline(1,'k-')
yline(1.5,'k-')

legend

%%
figure(1)
box off

xlabel('SEE length change (nm)')
ylabel('Force (-)')

%%
AMP = .0383 * newparms.s/2 * 1e9

