function[username, githubfolder] = get_paths()

usernames = {'timvd','u0167448'};

for i = 1:length(usernames)
    docfolder =  ['C:\Users\', usernames{i}, '\Documents'];
    
    if isfolder(docfolder)
        mainfolder = docfolder;
        username = usernames{i};
    end
end

if contains(mainfolder, 'timvd')
    githubfolder = mainfolder;
else
    githubfolder = [mainfolder, '\GitHub'];
end

addpath(genpath([githubfolder, '\muscle-thixotropy']))
addpath(genpath([mainfolder, '\casadi-3.7.1-windows64-matlab2018b']))
addpath(genpath([githubfolder, '\biophysical-muscle-model']))

end