clear all; close all; clc

[username, githubfolder] = get_paths();
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

iFs = 2;
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

%%   
clc
close all

X0 = parms.x0';

parms.vts = vis;
parms.ti = tis;

pCa     = 9;
Ca      = 10.^(-pCa+6);
parms.Cas = Ca;

ozs = [0 1];

figure(1)
color = get(gca,'colororder');

for j = 1:2
    parms.PE_isw_SE = ozs(j);

    sol = ode15s(@(t,y) fiber_dynamics_explicit_no_tendon(t,y, parms), [0 Ts(end)], X0, []);

    F = sol.y(1,:) + sol.y(2,:);
    Lce = sol.y(end-3,:);
    t = sol.x;
    
    if j == 1 % PE in parallel with SE
        L = interp1(tis, Lis * parms.gamma, t);
    else
        L = Lce; 
    end
    
    Fp = parms.Fpe_func(L, parms);
    Ftot = F * parms.Fscale + Fp;
    
    figure(1)
    subplot(211)
    plot(t, L, 'color', color(j,:)); hold on    
    
    subplot(212)
    plot(t, Ftot, 'color', color(j,:)); hold on
    plot(t, Fp, '--', 'color', color(j,:)); hold on

end

%%
