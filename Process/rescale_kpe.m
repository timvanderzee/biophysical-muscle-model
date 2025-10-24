clear all; close all; clc

mcode = [2 1 1];

if sum(mcode == [1 1 1]) == 3
    optparms = {'f', 'k11', 'k22', 'k21', 'J2', 'kon', 'kse', 'kse0'};
elseif sum(mcode == [1 1 3]) == 3
    optparms = {'f', 'k11', 'k22', 'k21', 'n','kappa', 'kse', 'kse0'};
elseif sum(mcode == [2 1 1]) == 3
    optparms = {'n','kappa', 'kse', 'kse0', 'vmax'};
%     optparms = {'n'};
end

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

casadifolder = 'C:\Users\u0167448\Documents\casadi-3.7.1-windows64-matlab2018b';

addpath(genpath([githubfolder, '\muscle-thixotropy\new_model\']))
addpath(genpath(casadifolder))
% rmpath(genpath('C:\Users\u0167448\Documents\casadi-3.7.1-windows64-matlab2018b'))
% addpath(genpath('C:\GBW_MyPrograms\casadi-3.6.7-windows64-matlab2018b'))
addpath(genpath([githubfolder, '\biophysical-muscle-model']))

%% specify data
load('active_trials.mat', 'Fm')
% iFs = [1 2 3, 5, 6, 7, 8, 10, 11];
iFs = [2,3,5,6,7,8,11];
iFs = 3;

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

for iF = iFs

% iF = 6; % fiber number (note: #4 and #9 lack pCa = 4.5)
Ks = find(Fm(:,iF) > .05); % only consider active trials
% Ks = find(isfinite(Fm(:,iF))); % only consider active trials
n = [3 1]; % ISI number
m = [7 1]; % AMP number
tiso = 3; % isometric time (s)

%% load data
cd(['C:\Users\',username,'\OneDrive - KU Leuven\9. Short-range stiffness\matlab\data'])
load([fibers{iF},'_cor_new.mat'],'data')

Data = prep_data_v2(data,n, m,Ks,tiso);

end

%%
t = data.texp(:,end,n(1),m(1));
F = data.Fexp(:,end,n(1),m(1));
L = data.Lexp(:,end,n(1),m(1));

%%
vs = {'\', '\'};
% vs = {'\full\'};

cd([githubfolder, '\muscle-thixotropy\new_model\get_variable'])
[output_mainfolder, filename, opt_type, ~] = get_folder_and_model(mcode);

disp(filename)
output_folder = [opt_type,'\normalized\with_PE_optimized\2_trials'];

output_dir = [output_mainfolder{1}, '\', filename,vs{1}, output_folder];
cd(output_dir)

% get other parameters from other fiber
% ii = 6; % fiber from which parameters are obtained
load([filename,'_F', num2str(iF),'_best.mat'],'parms','exitflag','fopt','C0','Cbounds','model','P0','P')
C = p_to_c(P, Cbounds);
parms = C_to_parms(C, parms, parms.optvars);
parms = calc_dependent_parms(parms);  

parms.act = 1;
parms.Noverlap = 1;
parms.approx = 0;
parms.K = 100;
parms.vF_func = @(vcerel,parms)parms.e(1)*log((parms.e(2)*vcerel./parms.vmax+parms.e(3))+sqrt((parms.e(2)*vcerel./parms.vmax+parms.e(3)).^2+1))+parms.e(4);

parms.gamma = .5*parms.s / parms.h; % length scaling
oparms = parms;

figure(1)
plot(t, F); hold on
plot(t, parms.Fpe_func(L * parms.gamma, parms))

%% adjust parameter
h = 10e-9; % powerstroke size
s = 2.6e-6; % sarcomere length

nparms = parms;
nparms.gamma = .5 * s ./ h;

nparms.kpe = parms.gamma / nparms.gamma * parms.kpe;

figure(1)
plot(t, parms.Fpe_func(L * nparms.gamma, nparms), '--')










