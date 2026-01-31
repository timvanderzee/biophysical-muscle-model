clear all; close all; clc

% this datafolder should contain the subfolders 2017 and 2018
datafolder = '';

% this folder should contain model fits
modelfolder = '';

% this folder should contain the biophysical-muscle-model repository
githubfolder = '';

%% Figures 2-4
for fig = 2:4
    cd([githubfolder, 'biophysical-muscle-model\Figures'])
    Figs2to4(datafolder, modelfolder, githubfolder, fig)
end

%% Figure 5
cd([githubfolder, 'biophysical-muscle-model\Figures'])
Fig5(githubfolder)

%% Figure 6
cd([githubfolder, 'biophysical-muscle-model\Figures'])
Fig6(githubfolder)

%% Figure 7
cd([githubfolder, 'biophysical-muscle-model\Figures'])
Fig7(datafolder, modelfolder, githubfolder)

%% Figure 8
cd([githubfolder, 'biophysical-muscle-model\Figures'])
Fig8(githubfolder)