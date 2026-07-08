function[datafolder, outputfolder] = get_paths(githubfolder)
    
%% output folder
% outputfolder is assumed to be in the same folder that the model code is
% in. If this folder does not exist, it is created. It will initially be
% empty. You'll have to run 'reproduce_model_forces' to fill this folder
% with model simulation output
outputfolder = fullfile(githubfolder, 'biophysical-muscle-model-results');

if ~isfolder(outputfolder)
    mkdir(outputfolder)
end

%% data folder
% datafolder is assumed to be in the same folder that the model code is in.
% If this folder does not exist, it is created. Note: at this time, data is
% provided to reviewers only. Reviewers will see the data as supplemental
% information to the submitted manuscript. Two subfolders are provided as
% ZIP files, named '2017' and '2018'. Please create a folder named
% 'biophysical-muscle-model-data' and put both unzipped folders (i.e. 2017
% and 2018) in this folder. If you are not a reviewer and do not have the
% data, you do not have to change the line below
datafolder = fullfile(githubfolder, 'biophysical-muscle-model-data');

if ~isfolder(datafolder)
    mkdir(datafolder)
end

end