clear all; close all; clc
cd(fileparts(which('reproduce_model_fitting.m')));
cd .. 
cd ..
githubfolder = cd;

% get folders
cd('biophysical-muscle-model')
addpath(genpath('Common'));
addpath(genpath('Reproduce'));
[datafolder, modelfolder] = get_paths(githubfolder);
save_results = 0;
visualize = 1;

%% step 1: specify which model we want to fit
model = '3-state XB coop'; 
[modelfunc, odefunc, modelname] = look_up_model(model);

%% step 2: specify which parameters we are fitting, and which bounds we are using
[optparms, bnds] = get_fitting_parms(model);

%% step 3: obtain (initial) parameter values
iF = 2; % [2,3,5,6,7,8,11];
[parms] = get_initial_parameters(model, bnds, githubfolder, iF);

%% step 4: specify which data we want to fit on
iF = 7; % [2,3,5,6,7,8,11];
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
N = 500; % number of interpolation points
th = .05; % threshold for active trials

if ishandle(1), close(1); end; figure(1)
[Xdata] = get_fitting_data(githubfolder, datafolder, iF, N, th, visualize);

%% step 5: obtain initial guess for model states
cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Fitting'))
IG = get_initial_guess(model, odefunc, modelfunc, Xdata, parms, visualize);

%% step 6: do fitting
% define weigth vector
w1 = 100;    % weight for fitting active force trajectory
w2 = 1000;   % weight for fitting passive force trajectory
w3 = 10; 	 % weight for regularization
weights = [w1 w2 w3];

% fit model parameters
[newparms, out, opti] = fit_model_parameters_v2(model, parms, optparms, bnds, Xdata, IG, weights);

%% step 7: analyze output
if ishandle(2), close(2); end
if ishandle(3), close(3); end
if ishandle(4), close(4); end

% load old parameters for reference
oldparms = load(fullfile(githubfolder, 'biophysical-muscle-model', 'Reproduce', 'Parameters',fibers{iF}, ['parms_', modelname, '.mat']), 'newparms');
if ~isfield(oldparms.newparms, 'Lce0')
    oldparms.newparms.Lce0 = -10;
end

figure(2)
for i = 1:length(optparms)
    nexttile
    bar(categorical({'IG', 'new', 'old'}), [parms.(optparms{i}) newparms.(optparms{i}) oldparms.newparms.(optparms{i})])
    title(optparms{i})
    box off
end

N = 1000;
th = .05;

figure(3)
[Xdata, Data] = get_fitting_data(githubfolder, datafolder, iF, N, th, visualize);
plot(out.t, out.F, '--')

for j = 1:2
    if j == 1
        testparms = newparms;
    else
        testparms = oldparms.newparms;
    end
    
    % use approximated exponentials
    testparms.approx = 1;

    % simulate
    [nsol, sout] = simulate_model(model, odefunc, modelfunc, Xdata, testparms);
    
    figure(3)
    subplot(414);
    plot(sout.t, sout.F,'-'); hold on
    
    % calc SRS
    [SRSrel] = calc_SRSrel(sout.t, sout.F, Xdata, Data, testparms);
    
    figure(4)
    plot(SRSrel.F0, SRSrel.ds(:,1)./SRSrel.ds(:,2),'o-'); hold on
    plot(SRSrel.F0, SRSrel.ns(:,1)./SRSrel.ns(:,2),'o-')
   
end

figure(3)
legend('Data', 'Input','Fit', 'Sim 1', 'Sim 2')

figure(4)
legend('Data', 'Sim 1', 'Data', 'Sim 2', 'location', 'best')
legend boxoff
box off
xlabel('Isometric force (F_0)')
ylabel('Relative stiffness')
xlim([0 1.05])


%% functions

function [Xdata, Data] = get_fitting_data(githubfolder, datafolder, iF, N, th, visualize)
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
n = [3 1]; % ISI number
m = [7 1]; % AMP number
tiso = 3; % isometric time (s)

% specify data
if isfolder(fullfile(datafolder, '2017')) % check whether there is a subfolder
    
    % load data
    years = {'2017', '2018'};
    for i = 1:length(years)
        if contains(fibers{iF}, years{i})
            fullfolder = fullfile(datafolder, years{i});
        end
    end
else
    fullfolder = datafolder; % otherwise, use the current folder
end


cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Data'))
load('active_trials.mat', 'Fm')

cd(fullfolder)
load([fibers{iF},'.mat'],'data')

Ks = find(Fm(:,iF) > th); % only consider active trials

cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Data', 'Processing'))
AData = prep_data_v2(data,n,m,Ks,tiso);

% passive
PData = prep_data_v2(data,3,7,length(data.pCas),tiso);

% combine
Data.t = [AData.t; PData.t(2:end) + AData.t(end)];
Data.F = [AData.F; PData.F(2:end)];
Data.L = [AData.L; PData.L(2:end)];
Data.v = [AData.v; PData.v(2:end)];
Data.C = [AData.C; PData.C(2:end)];

if visualize
    subplot(411)
    plot(Data.t, Data.C,'k.');
    
    subplot(412)
    plot(Data.t, Data.v,'k.');
    
    subplot(413)
    plot(Data.t, Data.L,'k.');
    
    subplot(414)
    plot(Data.t, Data.F,'k.');
    
    for j = 1:4
        subplot(4,1,j)
        box off
        hold on
    end
end

% convert to model input
[Atis, ACas, ALis, Avis, Ats] = create_input(tiso, AData.dTt, AData.dTc, AData.ISI, AData.Ca(Ks), N);
[Ptis, PCas, PLis, Pvis, Pts] = create_input(tiso, PData.dTt, PData.dTc, PData.ISI, PData.Ca(end), N);

% combine model input
tis = [Atis (Ptis(2:end) + Atis(end))];
Cas = [ACas PCas(2:end)];
Lis = [ALis PLis(2:end)];
vis = [Avis Pvis(2:end)];
ts = [Ats Pts + Ats(end)];

% interpolate force
Fis = interp1(Data.t, Data.F, tis);
Lis(Lis<0) = 0;

% intervals of interest
id1 = nan(length(Ks), length(tis), length(m));
id2 = nan(length(Ks), length(tis), length(m));
tmax = [0 tiso * length(Ks) tiso * (length(Ks)*2)];

for i = 1:length(m)
    for k = 1:length(Ks)
        id1(k,:,i) = (tis > (tmax(i) + ts(k) + tiso - 4*AData.dTt - 2*AData.dTc(i) - AData.ISI(i))) & (tis < (tmax(i) + ts(k)+tiso - AData.dTt));
        id2(k,:,i) = (tis > (tmax(i) + ts(k) + tiso - 5*AData.dTt - 2*AData.dTc(i) - AData.ISI(i))) & (tis < (tmax(i) + ts(k) + tiso - 4*AData.dTt - 2*AData.dTc(i) - AData.ISI(i)));
    end
end

% sum across conditions and pCas
idF = find(sum(sum(id1,3, 'omitnan'),1) & isfinite(Fis));
idC = find(sum(sum(id2,3, 'omitnan'),1));

% active indices (assume last trial is passive)
idA = 1:length(tis) - N;

% passive
idP = (idA(end)+1):length(tis);
idFP = find((tis > (tmax(3) + tiso - 4*PData.dTt - 2*PData.dTc - PData.ISI)) & (tis < (tmax(3) + tiso - PData.dTt)) & isfinite(Fis));

if visualize

    subplot(411)
    plot(tis, Cas, 'b');
    
    subplot(412)
    plot(tis, vis, 'b');
    
    subplot(413)
    plot(tis, Lis,'b');
    
    subplot(414)
    plot(tis, Fis, 'b');
    
end

% select data for fitting
Xdata.t = tis;
Xdata.F = Fis;
Xdata.v = vis;
Xdata.L = Lis;
Xdata.Cas = Cas;
Xdata.idF = idF(idF<idA(end));
Xdata.idC = [1:300 idC(idC<idA(end))];
Xdata.Ks = Ks;

% active

% if contains(model, 'XB')
Xdata.idA = idA;

% passive
Xdata.idP = idP;
Xdata.idFP = idFP;
% end

[id0,id1,id2] = get_indices(Data.t, tiso, ts, AData.dTt, AData.dTc(1), AData.ISI(1), AData.Ca(Ks));

Xdata.id0 = id0;
Xdata.id1 = id1;
Xdata.id2 = id2;

end

function[optparms, bnds] = get_fitting_parms(model)

% bounds
bnds.f = [20 1e3];
bnds.k11 = [1e-5 200];
bnds.k22 = [0 2];
bnds.k21 = [1 200];
bnds.kF = [1 1e4];
bnds.J1 = [1e-3 200];
bnds.J2 = [1 1e3];
bnds.kon = [5 200];
bnds.kse = [1e-3 1];
% bnds.kse0 = [1e-4 1];
bnds.kse0 = [1e-2 1];
bnds.koop = [1 200];
bnds.n = [.1 10];
bnds.kappa = [.1 10];
bnds.vmax = [1 200];
bnds.ps2 = [-1 2];
bnds.k = [1 5000];
bnds.b = [1 5000];
bnds.dLcrit = [1 10];
bnds.Lce0 = [-10 0];

bnds.kpe = [1e-4 1e-1];
bnds.Fpe0 = [1e-5 1e-1];

% parameters to be fitted
if strcmp(model, '3-state XB coop')
    optparms = {'f', 'k11', 'k22', 'k21', 'kon', 'kse', 'kse0', 'kF', 'koop', 'Lce0', 'kpe'};
    
elseif strcmp(model, '2-state XB coop')
    optparms = {'f', 'k11', 'k22', 'k21', 'kon', 'kse', 'kse0', 'koop','Lce0', 'kpe'};
    
elseif strcmp(model, '2-state XB')
    optparms = {'f', 'k11', 'k22', 'k21', 'n','kappa', 'kse', 'kse0','Lce0', 'kpe'};
    
elseif strcmp(model, 'Hill-type SE')
    optparms = {'n','kappa', 'kse', 'kse0', 'vmax', 'Lce0', 'kpe'};
    
elseif strcmp(model, 'Hill-type no SE')
    optparms = {'n','kappa', 'vmax', 'Lce0', 'kpe'};
    
elseif strcmp(model, '4-state XB coop')
    optparms = {'f', 'k11', 'k22', 'k21', 'kon', 'kse', 'kse0', 'kF', 'koop', 'dLcrit', 'Lce0', 'kpe'};
end


end

function[parms] = get_initial_parameters(model, bnds, githubfolder, iF)

[~, ~, modelname] = look_up_model(model);

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Reproduce', 'Process'))
% [~, ~, modelname] = get_model_folder('', mcode, 0);

foldername = fullfile(githubfolder, 'biophysical-muscle-model', 'Reproduce', 'Parameters',fibers{iF});

cd(foldername)
load(['parms_', modelname, '.mat'], 'redparms')
newparms = complete_parms(redparms);
parms = newparms;

% add some parameters that didn't exist before
parms.PE_isw_SE = 1;
parms.Fse_func = @(dlse, parms) parms.kse0*(exp(parms.kse*dlse)-1);
parms.Fpe0 = parms.Fpe0/2;
parms.Lce0 = -1;
parms.Fpe_func = @(L, parms) parms.kpe*(L-parms.lmtc0).*(L>parms.lmtc0)+parms.Fpe0;
parms.kpe_func = @(Lce, parms) parms.kpe;

parms.kF = min(parms.kF, .9*bnds.kF(2));
parms.kse0 = max(parms.kse0, 1.1*bnds.kse0(1));

if ~isfield(parms, 'gamma')
    parms.gamma = .5*parms.s / parms.h; % length scaling
end

% store old parameters
oparms = parms;

% update powerstroke
h = 10e-9; % powerstroke size
s = 2.6e-6; % sarcomere length
parms.gamma = .5*s / h; % length scaling
parms.kpe = oparms.gamma / parms.gamma * oparms.kpe;

% approximate
parms.approx = 0;

if strcmp(model, '3-state XB coop')
    %         parms.J1 = 6.17;
    %         parms.koop = 5.7;
    %         parms.JF = 1e3;
    %         parms.J2 = 200;
    
elseif strcmp(model, '2-state XB coop')
    parms.J1 = 0;
    parms.J2 = 0;
    parms.JF = 0;
    parms.kF = 0;
    
elseif strcmp(model, 'Hill-type SE')
    parms.f = 0;
    parms.k = 0;
    parms.dLcrit = 0;
    
    parms.act_max = (1 - parms.Fpe0) / parms.Fscale;
    
elseif strcmp(model, '4-state XB coop')
    %         parms.J1 = 6.17;
    %         parms.koop = 5.7;
    %         parms.JF = 1e3;
    %         parms.J2 = 200;
    
    parms.k = 1000;
    parms.b = 3000;
    
    %         parms.k = 500;
    %         parms.b = 1000;
    parms.dLcrit = 2;
    
    parms.ps2 = 0;
end
end

function[SRSrel, RMSD] = calc_SRSrel(t, Fse, Xdata, Data, parms)
    
% interpolate
nFi = interp1(t, Fse, Data.t);
Lti = interp1(Xdata.t, Xdata.L * parms.gamma, Data.t);

% compute RMSD
RMSD = mean(abs(nFi - Data.F));

% pre-allocate
SRSrel.ns = nan(length(Xdata.Ks), 2);
SRSrel.ds = nan(length(Xdata.Ks), 2);
SRSrel.F0 = nan(length(Xdata.Ks), 1);

for i = 1:length(Xdata.Ks)

    np1 = polyfit(Lti(Xdata.id1(i,:)), nFi(Xdata.id1(i,:)), 1);
    np2 = polyfit(Lti(Xdata.id2(i,:)), nFi(Xdata.id2(i,:)), 1);

    dp1 = polyfit(Data.L(Xdata.id1(i,:)), Data.F(Xdata.id1(i,:)), 1);
    dp2 = polyfit(Data.L(Xdata.id2(i,:)), Data.F(Xdata.id2(i,:)), 1);

    SRSrel.ds(i,:) = [dp1(1) dp2(1)];
    SRSrel.ns(i,:) = [np1(1) np2(1)];
    SRSrel.F0(i) = mean(Data.F(Xdata.id0(i,:)));

end
end