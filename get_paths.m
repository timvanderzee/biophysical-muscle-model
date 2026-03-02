function[datafolder, outputfolder] = get_paths(githubfolder)
    
%% output folder
% datafolder should be a string specifying the location of the folder where
% you would like to export the simulation output to. It can be any folder
% on your device.
outputfolder = fullfile(githubfolder, 'biophysical-muscle-model-results');

if ~isfolder(outputfolder)
    mkdir(outputfolder)
end

%% data folder
% datafolder should be a string specifying the location of the folder
% containing the experimental data. Note: at this time, data is provided to
% reviewers only. Reviewers will see the data as supplemental information
% to the submitted manuscript. Two subfolders are provided as ZIP files,
% named '2017' and '2018'. datafolder should specify the folder that
% contains both these (unzipped) subfolders. If you are not a reviewer and
% do not have the data, you do not have to change the line below
datafolder = fullfile(githubfolder, 'biophysical-muscle-model-data');

if ~isfolder(datafolder)
    mkdir(datafolder)
end

end