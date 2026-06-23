clear all; close all; clc
cd(fileparts(which('reproduce_manuscript_figures.m')));
cd .. 
cd ..
githubfolder = cd;

% get folders
cd('biophysical-muscle-model')
[datafolder, modelfolder] = get_paths(githubfolder);

%% Figures 2-4: Force traces
for fig = 2:4
    cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Reproduce', 'Figures'))
    Figs2to4(datafolder, modelfolder, githubfolder, fig)
end

%% Figure 5: Force RMSD
cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Reproduce', 'Figures'))
Fig5(githubfolder)

%% Figure 6: SRS
cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Reproduce', 'Figures'))
Fig6(githubfolder)

%% Figure 7: SRS
cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Reproduce', 'Figures'))
Fig7(datafolder, modelfolder, githubfolder)

%% Figure 8: Computational cost
cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Reproduce', 'Figures'))
Fig8(githubfolder)