clear all; close all; clc

%% main data folder
% containing the 'raw' data files
% note: these are currently not provided
datafolder = '';
current_folder = cd; % needs to be the folder where this script lives
addpath(genpath(fullfile(current_folder, 'Processing')))

%% reorganize stretch-shortening data
raw_datafolder = fullfile(datafolder, 'Transformed Raw Data to Matlab');
reorganized_datafolder = fullfile(datafolder, 'Reorganized Data');

reorganize_data(raw_datafolder, reorganized_datafolder)
% result: e.g., 6Aug2018a_reorganized.mat

%% get force-pCa from activation trials
activation_datafolder = fullfile(datafolder, 'Activation Data');
force_pCa_folder = current_folder;

calc_force_pCa(activation_datafolder, force_pCa_folder)
% result: force_pCa.mat

%% normalize stretch-shortening data based on force-pCa
normlized_datafolder = fullfile(datafolder, 'Normalized Data');
normalize_data(force_pCa_folder, reorganized_datafolder, normlized_datafolder)
% result: e.g., 6Aug2018a.mat