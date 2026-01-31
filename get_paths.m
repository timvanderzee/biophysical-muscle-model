function[datafolder, modelfolder, githubfolder] = get_paths()
    
% this datafolder should contain the subfolders 2017 and 2018
datafolder = '';

% this folder is where the model fits are or will be saved
modelfolder = '';

% this folder should contain the biophysical-muscle-model repository
githubfolder = '';

% add the model functions to the path
addpath(genpath([githubfolder, '\biophysical-muscle-model\Common\']))

end