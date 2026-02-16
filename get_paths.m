function[datafolder, outputfolder, githubfolder] = get_paths()
    
%% github folder
% githubfolder should be a string specifying the location of the folder
% containing the biophysical-muscle-model repository. It depends on where
% you cloned or downloaded the repository to
githubfolder = 'C:\Users\u0167448\Documents\GitHub'; % e.g. C:\John_Doe\Documents\github

% add the model functions to the path
addpath(genpath([githubfolder, '\biophysical-muscle-model\Common\']))

%% output folder
% datafolder should be a string specifying the location of the folder where
% you would like to export the simulation output to. It can be any folder
% on your device.
outputfolder = 'C:\Users\u0167448\Desktop\Model fits'; % e.g. C:\John_Doe\Documents\simulation_output

%% data folder
% datafolder should be a string specifying the location of the folder
% containing the experimental data. Note: at this time, data is provided to
% reviewers only. Reviewers will see the data as supplemental information
% to the submitted manuscript. Two subfolders are provided as ZIP files,
% named '2017' and '2018'. datafolder should specify the folder that
% contains both these (unzipped) subfolders. If you are not a reviewer and
% do not have the data, you do not have to change the line below
datafolder = 'C:\Users\u0167448\OneDrive - KU Leuven\Horslen paper\Normalized Data'; % e.g. C:\John_Doe\Documents\data

end