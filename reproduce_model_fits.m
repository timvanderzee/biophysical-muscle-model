clear all; close all; clc

% get folders
[~, modelfolder, githubfolder] = get_paths();

%% Choose the fiber
% valid options: (2,3,5,6,7,8,11);
% note: 1, 4, 9 and 10 were excluded (see pre-print)
fiber = 6;

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

%% Calculate model time-series
% note: you could insert nested for loops across fibers, models and methods
% here to obtain all model predictions
cd([githubfolder, '\biophysical-muscle-model\Process'])
visualize = 0;
calc_model_forces(githubfolder, modelfolder, model, fiber, discretized, visualize)

