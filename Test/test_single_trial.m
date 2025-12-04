clear all; close all ;clc
[username, githubfolder] = get_paths();
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

iFs = 6;
parms_version = '_v3';

% load parameters
mcode = [1 1 1];
[output_mainfolder, modelname, ~, ~] = get_folder_and_model(mcode);

input_foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iFs}];
cd(input_foldername)
load(['parms_',modelname, parms_version, '.mat'], 'newparms')
parms = update_parms(newparms);


pCa     = 9;
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
parms.Lts = [Lis Lis Lis Lis] * parms.gamma;

parms.x0 = zeros(1,6);

parms.Fpe_func = @(Lce, parms) parms.kpe * log(1+exp(Lce*parms.K))/parms.K + parms.Fpe0;
parms.kpe_func = @(Lce, parms) parms.kpe .* (1 - 1./(exp(parms.K*Lce)+1));

% parms.Fpe_func = @(L, parms) parms.kpe*(L-parms.lmtc0).*(L>parms.lmtc0)+parms.Fpe0;
% parms.kpe_func = @(Lce, parms) parms.kpe;

Q0 = 0;
p0 = 0;
q0 = .1;


%%
close all
odeopt = odeset('maxstep', 1e-1);


parms.K = 1e2;

for j = [1 3]
    if j == 1
        [Q00, Q20, lce0, Q10] = find_steady_state(Q0, p0, q0, parms, 'regular');
        parms.x0 = [Q00 Q20 lce0 0 0 0];
        sol = ode15s(@(t,y) fiber_dynamics_explicit_length(t,y, parms), [0 max(parms.ti)], parms.x0(:), odeopt);
    elseif j == 2
        [Q00, Q20, lce0, Q10] = find_steady_state(Q0, p0, q0, parms, 'adjusted');
        parms.x0 = [Q00 Q10 Q20 lce0 0 0 0];
        sol = ode15s(@(t,y) fiber_dynamics_explicit_no_tendon(t,y, parms), [0 max(parms.ti)], parms.x0(:), odeopt);
    elseif j == 3
        [Q00, Q20, lce0, Q10] = find_steady_state(Q0, p0, q0, parms, 'regular');
        parms.x0 = [Q00 Q20 lce0 0 0 0];
        parms.x0p = zeros(size(parms.x0));
        sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_length(t,y,yp, parms), [0 max(parms.ti)], parms.x0(:), parms.x0p(:), odeopt);
    end
    
    
    L = sol.y(end-3,:);
    t = rem(sol.x, max(Ts)); 
    t = sol.x;
    dlse = interp1(parms.ti, parms.Lts, sol.x) - L;
    
    parms.Fse_func = @(dlse, parms) parms.kse0*(exp(parms.kse*dlse)-1);
    Fse = parms.Fse_func(dlse, parms);
     
    if j == 2
        Fpe = parms.Fpe_func(interp1(parms.ti, parms.Lts, sol.x), parms);
        F = Fse + Fpe;
    else
        F = Fse;
    end

    figure(1)

    subplot(311)
    plot(t, L+dlse); hold on

    subplot(312)
    plot(t, L); hold on

    subplot(313)
    plot(t, F); hold on
end

    % plot(t, F1-F2)

%%
N = length(sol.x);
id = N;

dy = fiber_dynamics_explicit_length(sol.x(id), sol.y(:,id), parms)



