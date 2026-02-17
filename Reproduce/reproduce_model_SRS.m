clear all; close all; clc
cd ..

% get folders
[~, modelfolder, githubfolder] = get_paths();

%% Choose the model
% 1. Hill-type without SEE
% 2. Hill-type with SEE
% 3. Biophysical without cooperative dynamics
% 4. Biophysical with only thin filament cooperative dynamics
% 5. Biophysical with both thin and thick filament cooperative dynamics
% 6. Biophysical with cooperative dynamics and forcibly detached state
model = 5;

%% Choose discretized or approximated solution method
% 1. Discretized solution method (method of characteristics)
% 2. Approximated solution method (distribution-moment)
% note: only applicable if a biophysical model has been selected
discretized = 0;

%% Reproduce model SRS
cd(fullfile(githubfolder, 'biophysical-muscle-model','Reproduce', 'Process'))
calc_model_SRS(githubfolder, modelfolder, model, discretized)