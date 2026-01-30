clear all; close all; clc
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
[username, githubfolder] = get_paths();
iFs = [2,3,5,6,7,8,11];

versions = {'', '', '_v2', '', '_v4', '_v4'};

mcodes = [2 2 1; 2 1 1; 1 1 3; 1 1 2; 1 1 1; 1 2 1];
optparms = {'n','kappa', 'vmax', 'f', 'k11', 'k22', 'k21', 'kon', 'kse', 'kse0', 'JF', 'koop', 'Lce0', 'kpe', 'dLcrit'};

for j = 1:size(mcodes,1)
    mcode = mcodes(j,:);
    [output_mainfolder, filename, opt_type, ~] = get_folder_and_model(mcode);
    disp(filename)

    for i = 1:length(iFs)

        foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iFs(i)}];

        cd(foldername)
        load(['parms_', filename, versions{j}, '.mat'], 'newparms')

        for k = 1:length(optparms)
            if isfield(newparms, optparms{k})
                eval([optparms{k}, '(j,i) = ', num2str(newparms.(optparms{k})),';'])
            else
                eval([optparms{k}, '(j,i) = NaN'])
            end
                
        end
        
        
        % evaluate parameters
%         f(j,i) = newparms.f;

    end
end

%% 
for k = 1:length(optparms)
   disp(optparms{k})
   eval(['disp(median(', optparms{k}, ',2))'])
end
