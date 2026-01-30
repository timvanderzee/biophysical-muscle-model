clear all; close all; clc

% this folder should contain model fits
modelfolder = 'C:\Users\u0167448\Desktop\Model fits';

% this folder should the biophysical-muscle-model repository
githubfolder = 'C:\Users\u0167448\Documents\GitHub\';

%% Calculate model time-series
calc_model_forces(githubfolder, modelfolder)

