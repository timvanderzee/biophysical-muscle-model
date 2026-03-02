clear all; close all; clc
cd(fileparts(which('create_reduced_parms.m')));

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

allparms = {'f', 'b', 'k11', 'k12', 'k21', 'k22', 'J1', 'J2', 'JF', 'kon', 'koff', 'koop', 'Noverlap', 'Fscale', 'w', 'Fpe0', 'approx', ...
    'kse', 'kse0', 'kpe', 'lce0', 'vmax', 'e', 'act_max', 'n', 'kappa', 'k', 'dLcrit', 'x0', 'amin', 'K', 'gamma', 'kF', 'Lce0', 'xi'};


cd(fibers{2})

files = dir(cd);

for i = 1:length(files)-2
    redparms = struct();

    clearvars bnds newparms optparms out
    load(files(i+2).name, 'bnds', 'newparms', 'optparms', 'out')
    
    newparms.xi = linspace(-15,15,500);
    
    for j = 1:length(allparms)
        if isfield(newparms, allparms{j})
            redparms.(allparms{j}) = newparms.(allparms{j});
        end
    end
    
    save(files(i+2).name, 'redparms', '-append')
end
    