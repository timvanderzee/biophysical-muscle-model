clear all; close all ;clc
[username, githubfolder] = get_paths();
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

iFs = 7;
parms_version = '_v3';

% load parameters
mcode = [1 1 1];
[output_mainfolder, modelname, ~, ~] = get_folder_and_model(mcode);

input_foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iFs}];
cd(input_foldername)
load(['parms_',modelname, parms_version, '.mat'], 'newparms')
parms = update_parms(newparms);


pCa     = 6.1;
Ca      = 10.^(-pCa+6);
AMP     = .0383;
ISI     = .001;
dTt     = .0383/.4545; % test stretch (= constant)
dTc     = AMP / .4545; % conditioning stretch
tiso    = dTt*3+dTc*2+ISI + 2;
dt      = .001; % gives 10 points in SRS zone
N       = round(tiso / dt);
[tis, Cas, Lis, vis, ts, Ts] = create_input(tiso, dTt, dTc, ISI, Ca, N);

% load('protocol.mat')

%%
clc
parms.vts = [vis vis vis vis];
parms.ti = linspace(0, max(tis)*4, length(parms.vts));
parms.Cas = [Cas Cas*.01 Cas Cas];
parms.Lis = [Lis Lis Lis Lis] * parms.gamma;

parms.x0(1:3) = .05;
F0 = parms.x0(1)+parms.x0(2) + parms.Fpe0;
dlse0 = parms.Lse_func(F0, parms);
lce0 = -dlse0;
parms.x0(end-3) = lce0;

% parms.vts = vts;
% parms.Lts = Lts;
% parms.Cas = Cas;
% parms.ti = toc;

odeopt = odeset('maxstep', 1e-2);

parms.K = 5e2;
sol = ode15s(@(t,y) fiber_dynamics_explicit_no_tendon(t,y, parms), [0 max(parms.ti)], [parms.x0(:); 0], odeopt);

%%
close all
L = sol.y(end-4,:);
t = sol.x; 
% Fp = parms.Fpe_func(L, parms);

figure(1)
Fce = sol.y(1,:) + sol.y(2,:);
% F1 = Fce;
% k   = parms.K;
F1   = Fce + parms.Fpe_func(L, parms);

% dlse = interp1(parms.ti, parms.Lis, sol.x) - L;
dlse = sol.y(end,:)*parms.gamma - L;

parms.Fse_func = @(dlse, parms) parms.kse0*(exp(parms.kse*dlse)-1);
F2 = parms.Fse_func(dlse, parms);

subplot(311)
plot(t, sol.y(end,:)*parms.gamma, '-', t, L+dlse, '--')

subplot(312)
plot(t, L, t, dlse)

subplot(313)
plot(t, [F1; F2])
% plot(t, F1-F2)




