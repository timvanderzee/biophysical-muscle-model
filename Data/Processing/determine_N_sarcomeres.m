clear all; close all; clc

cd('C:\Users\u0167448\OneDrive - KU Leuven\Horslen paper\Transformed Raw Data to Matlab')

files = dir('*.mat');

fL = nan(100,1);
sL = nan(100,1);

for i = 1:length(files)
    clear new_data
    load(files(i).name, 'new_data');
    
    if exist('new_data', 'var')
        sL(i) = new_data(1).sarcomere_length;
        fL(i) = new_data(1).muscle_length;
    end
end

%%
sLs = sL(:);

% assumed a typo
sLs(sLs > 1e-5) = sLs(sLs > 1e-5)/10;

figure(10)
histogram(sLs)

N = fL(:) ./ sLs;