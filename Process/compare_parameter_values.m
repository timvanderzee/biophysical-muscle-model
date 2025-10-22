clear all; close all; clc
[username, githubfolder] = get_paths();

mcodes = [2 1 1; 1 1 1; 1 1 3; 1 2 1];
mcodes = [1 1 1];

% iFs = [1 2 3, 5, 7, 8, 10, 11];
% iFs = 6;
iFs = [2,3,5,6,7,8,11];

AMPs = [0 12 38 121 216 288 383 682]/10000;
% ISIs = [1 10 100 316 1000 3160 10000]/1000;
ISIs = [.05 .2 .5];
% AMPs = [0    0.0012    0.0038    0.0121    0.0216    0.0288    0.0383    0.0532    0.0682];
% ISIs = [ 0.0010    0.0100    0.0500    0.1000    0.2000    0.3160    0.5000    1.0000    3.1600   10.0000];

pCas = [4.5 6.1 6.2 6.3 6.4 6.6 9];
Ca = 10.^(-pCas+6);
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

visualize = 0;

version = '_v2';

% AMPs = 682 / 10000;
% ISIs = [3160 10000]/1000;
% 


% load parameters
mcode = [1 1 1];
[output_mainfolder, modelname, ~, ~] = get_folder_and_model(mcode);

%% evaluate
for iF = iFs
    
    % disp(filename)
    foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iF}];
    cd(foldername)
    load(['parms_',modelname, '.mat'], 'newparms')
    
    parms(iF) = newparms;

end