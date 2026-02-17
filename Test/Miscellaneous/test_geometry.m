clear all; close all; clc

[username, githubfolder] = get_paths();
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

iFs = 8

for iF = iFs
parms_version = '_v3';

% load parameters
mcode = [1 1 1];
[output_mainfolder, modelname, ~, ~] = get_folder_and_model(mcode);

input_foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iF}];
cd(input_foldername)
load(['parms_',modelname, parms_version, '.mat'], 'newparms')
parms = update_parms(newparms);

pCa     = 6.3;
Ca      = 10.^(-pCa+6);
AMP     = .0383;
ISI     = .001;
dTt     = .0383/.4545; % test stretch (= constant)
dTc     = AMP / .4545; % conditioning stretch
tiso    = dTt*3+dTc*2+ISI + 2;
dt      = .001; % gives 10 points in SRS zone
N       = round(tiso / dt);
[tis, Cas, Lis, vis, ts, Ts] = create_input(tiso, dTt, dTc, ISI, Ca, N);

X0 = parms.x0';

parms.vts = vis;
parms.ti = tis;
parms.Cas = Cas;

ozs = [0 1 1];
ls = {'-','--',':'};

figure(1)
color = get(gca,'colororder');


for j = 1:3
    parms.PE_isw_SE = ozs(j);
    
    if j < 3 % using old derivative function
        sol = ode15s(@(t,y) fiber_dynamics_explicit_no_tendon(t,y, parms), [0 Ts(end)], X0, []);
        F = sol.y(1,:) + sol.y(2,:);
        
        Lce = sol.y(end-3,:);
        t = sol.x;

        if j == 1 % PE in parallel with SE
            Ltot = interp1(tis, Lis * parms.gamma, t);
         	Fp = parms.Fpe_func(Ltot, parms);
        else
            Fp = parms.Fpe_func(Lce, parms);
        end

        Ftot = F * parms.Fscale + Fp;
        
    else
    
%         [Q00, Q20, lce0, Q10] = find_steady_state(0, -1, 0, parms, 'regular');
        
        % find dlse where Fse = Fpe (assuming Fce = 0)
        dlse = parms.Lse_func(parms.Fpe0, parms);
        
        lce0 = -dlse;
        X0 = [1e-6 -1e-6 lce0 0 0 0 0];
        
        parms.lmtc0 = lce0; 
        
        sol = ode15s(@(t,y) fiber_dynamics_explicit_length(t,y, parms), [0 Ts(end)], X0, []);
        
        Lce = sol.y(3,:);
        Fp = parms.Fpe_func(Lce, parms);
        t = sol.x;
        Ltot = sol.y(end,:);
        dlse = Ltot - Lce;
        
        % note: -Fp needed because it was previously not scaled
        Ftot = parms.Fse_func(dlse, parms) * parms.Fscale - Fp;
    end
    
        figure(iF)
        subplot(311)
        plot(t, Lce, 'color', color(j,:), 'linestyle', ls{j}, 'linewidth',2); hold on
        
        subplot(312)
        plot(t, Ftot, 'color', color(j,:), 'linestyle', ls{j}, 'linewidth',2); hold on
        
        subplot(313);
        plot(t, Fp, 'color', color(j,:), 'linestyle', ls{j}, 'linewidth',2); hold on
        
end
    
end
