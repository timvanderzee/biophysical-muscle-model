clear all; close all; clc
[username, githubfolder] = get_paths();
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

% model to be fitted
mcode = [1 1 1];
version = '_v3';
% version = '';

% settings
N = 500;
% N = 1500;
save_results = 0;
visualize = 1;

iFs = 11 %; [2,3,5,6,7,8,11];
% n = [3 1]; % ISI number
% m = [7 1]; % AMP number
n = [3]; % ISI number
m = [7]; % AMP number
tiso = 3; % isometric time (s)

% bounds
bnds.Lce0 = [-100 0];
bnds.kpe = [1e-4 1e-1];
% parameters to be fitted

optparms = {'Lce0', 'kpe'};

%% specify data
load('active_trials.mat', 'Fm')

iF = iFs;

%% load data
cd(['C:\Users\',username,'\OneDrive - KU Leuven\9. Short-range stiffness\matlab\data'])
load([fibers{iF},'_cor_new.mat'],'data')

Ks = find(Fm(:,iF) < 0.05); % only consider active trials
Data = prep_data_v2(data,n, m,Ks,tiso);
[tis, Cas, Lis, vis, ts] = create_input(tiso, Data.dTt, Data.dTc, Data.ISI, Data.Ca(Ks), N);

Lis(Lis<0) = 0;

% interpolate force
Fis = interp1(Data.t, Data.F, tis);

if visualize
    if ishandle(1), close(1); end
    figure(1)
    subplot(411)
    plot(Data.t, Data.C,'k.', tis, Cas, 'b');
    
    subplot(412)
    plot(Data.t, Data.v,'k.', tis, vis, 'b');
    
    subplot(413)
    plot(Data.t, Data.L,'k.', tis, Lis,'b');
    
    subplot(414)
    plot(Data.t, Data.F,'k.');
    
    for j = 1:4
        subplot(4,1,j)
        box off
        hold on
    end
end

cd([githubfolder, '\muscle-thixotropy\new_model\get_variable'])
[output_mainfolder, filename, opt_type, ~] = get_folder_and_model([1 1 1]);

disp(filename)

foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iF}];
cd(foldername)
load(['parms_', filename, version, '.mat'], 'newparms')

parms = newparms;

%% initial guess: solve force equilibrium
parms.ti = tis;
parms.vts = vis;
parms.Cas = Cas;
parms.Lts = Lis * parms.gamma;

parms.Lce0 = 0;

Lcost = @(dlse, Lis, parms) ((parms.kse0*(exp(parms.kse*dlse)-1)) - (parms.kpe * (Lis*parms.gamma - dlse - parms.Lce0))).^2;

for i = 1:length(Lis)
    dlse(i) = fminsearch(@(L) Lcost(L, Lis(i), parms), 0);
end

dLce = (Lis * parms.gamma - dlse) - parms.Lce0;

Fse = parms.Fse_func(dlse, parms) * parms.Fscale;
Fpe = parms.kpe * dLce  * parms.Fscale;

figure(1)
subplot(414)
plot(tis, Fse, '-', tis, Fpe, '--')

% initial guess
IG.F = Fse;

%% intervals of interest
id1 = nan(length(Ks), length(parms.ti), length(m));
id2 = nan(length(Ks), length(parms.ti), length(m));
tmax = [0 tiso * length(Ks)];

for i = 1:length(m)
    for k = 1:length(Ks)
        id1(k,:,i) = (parms.ti > (tmax(i) + ts(k) + tiso - 4*Data.dTt - 2*Data.dTc(i) - Data.ISI(i))) & (parms.ti < (tmax(i) + ts(k)+tiso - Data.dTt));
        id2(k,:,i) = (parms.ti > (tmax(i) + ts(k) + tiso - 5*Data.dTt - 2*Data.dTc(i) - Data.ISI(i))) & (parms.ti < (tmax(i) + ts(k) + tiso - 4*Data.dTt - 2*Data.dTc(i) - Data.ISI(i)));
    end
end

% sum across conditions and pCas
idF = find(sum(sum(id1,3, 'omitnan'),1) & isfinite(Fis));
idC = find(sum(sum(id2,3, 'omitnan'),1));

if visualize
    subplot(414); hold on
    plot(tis(idF), Fse(idF),'m*', 'markersize', 1); hold on
    plot(tis(idC), Fse(idC),'g.'); hold on
end

%% select data for fitting
Xdata.t = tis;
Xdata.F = Fis;
Xdata.v = vis;
Xdata.L = Lis * parms.gamma;
Xdata.Cas = Cas;
Xdata.idF = idF;
Xdata.idC = idC;

%% fit!
[newparms, out, opti] = fit_model_parameters_passive(optparms, [], Xdata, parms, IG, bnds);

%% test
for i = 1:length(Lis)
    dlse(i) = fminsearch(@(L) Lcost(L, Lis(i), newparms), 0);
end

dLce = (Lis * newparms.gamma - dlse) - newparms.Lce0;

nFse = newparms.Fse_func(dlse, newparms) * newparms.Fscale;
nFpe = newparms.kpe * dLce  * newparms.Fscale;

% figure(1)
plot(Data.t, Data.F,'k.',tis, nFse, '-', tis, nFpe, '--', out.t, out.F, ':')









