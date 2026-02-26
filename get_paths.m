function[datafolder, outputfolder, githubfolder] = get_paths()
    
%% github folder
% githubfolder should be a string specifying the location of the folder
% containing the biophysical-muscle-model repository. It depends on where
% you cloned or downloaded the repository to
cd ..
githubfolder = cd; % e.g. C:\John_Doe\Documents\github

% add the model functions to the path
addpath(genpath([githubfolder, '\biophysical-muscle-model\Common\']))

if ~isfolder(githubfolder)
    warning('githubfolder does not point to a folder')
end

%% output folder
% datafolder should be a string specifying the location of the folder where
% you would like to export the simulation output to. It can be any folder
% on your device.
outputfolder = 'test'; % e.g. C:\John_Doe\Documents\simulation_output

if ~isfolder(outputfolder)
    warning('outputfolder does not point to a folder')
end

%% data folder
% datafolder should be a string specifying the location of the folder
% containing the experimental data. Note: at this time, data is provided to
% reviewers only. Reviewers will see the data as supplemental information
% to the submitted manuscript. Two subfolders are provided as ZIP files,
% named '2017' and '2018'. datafolder should specify the folder that
% contains both these (unzipped) subfolders. If you are not a reviewer and
% do not have the data, you do not have to change the line below
datafolder = 'test'; % e.g. C:\John_Doe\Documents\data

if ~isfolder(datafolder)
    warning('datafolder does not point to a folder')
end

end