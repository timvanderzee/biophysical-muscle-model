clear all; close all; clc
cd ..

% get folders
[datafolder, modelfolder, githubfolder] = get_paths();
save_results = 0;
visualize = 1;

%% step 1: specify which model we want to fit
model = '3-state XB coop'; 
[modelfunc, modelname] = look_up_model(model);

%% step 2: specify which parameters we are fitting, and which bounds we are using
[optparms, bnds] = get_fitting_parms(model);

%% step 3: obtain (initial) parameter values
iF = 2; % [2,3,5,6,7,8,11];
[parms] = get_initial_parameters(model, bnds, githubfolder, iF);

%% step 4: specify which data we want to fit on
% create a given ramp protocol
T = 1;          % total duration (s)
dt = 1/1000;    % sample time (s)
Ca = 1;         % calcium concentration (uM)
RT = .1;        % recovery time(s)
A = .04;        % stretch amplitude (L0)

Xdata = ramp(Ca, T, dt, RT, A);
color = lines(3);

% let's create some fake data to fit. in this example, we start with the
% model force corresponding to the above input (Xdata), the model (defined
% in step 1), and the model parameters (defined in step 3) we manipulate
% this force by multiply it with a constant, and adding an offset and noise
cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Fitting'))
IG = get_initial_guess(Xdata, parms, 0);
Xdata.F = .2 + IG.Fsei * 1.3 + randn(size(IG.Fsei)) * .01;

% visualize data to be fitted
if ishandle(1), close(1); end
figure(1)
subplot(411)
plot(Xdata.t, Xdata.Cas, 'color', color(1,:), 'linewidth', 1.5); hold on; box off
ylabel('[Ca]^{2+} (\muM)')
title('Calcium concentration')

subplot(412)
plot(Xdata.t, Xdata.L, 'color', color(1,:), 'linewidth', 1.5); hold on; box off
ylabel('\DeltaFiber length (L_0)')
title('Fiber length change')

subplot(413)
plot(Xdata.t, Xdata.v, 'color', color(1,:), 'linewidth', 1.5); hold on; box off
ylabel('Fiber velocity (L_0/s)')
title('Fiber velocity')

subplot(414)
plot(Xdata.t, Xdata.F, 'color', color(1,:), 'linewidth', 1.5); hold on; box off
ylabel('Fiber force (F_0)')
title('Fiber force')
    
%% step 5: obtain initial guess for model states
% we have already obtained the initial guess (IG) corresponding to the
% model, parameters and input in the previous step

subplot(414)
plot(Xdata.t, IG.Fsei, 'color', color(2,:), 'linewidth', 1.5); hold on; box off
legend('Target', 'Initial guess', 'location', 'best')

%% step 6: do fitting
% define weight vector
weights = 100;

% fit model parameters
[newparms, out, opti] = fit_model_parameters_v2(model, parms, optparms, bnds, Xdata, IG, weights);

%% step 7: analyze output
% get initial state
X0 = [out.Q0(1) out.Q2(1) out.Fse(1) out.Non(1) out.DRX(1) 0];

% simulate
newparms.ti = Xdata.t;
newparms.vts = Xdata.v;
newparms.Cas = Xdata.Cas;
newparms.Lts = Xdata.L * parms.gamma;

nsol = ode15s(@(t,y,yp) modelfunc(t,y, newparms), [0 max(newparms.ti)], X0, odeset('maxstep', 1e-2));

% get force
[t, Fse, Fpe, Fce] = get_forces_from_state(nsol, Xdata, newparms);  

%% visualize output
if ishandle(2), close(2); end
figure(2)
for i = 1:length(optparms)
    nexttile
    bar(categorical({'old', 'new'}), [parms.(optparms{i}) newparms.(optparms{i})])
    title(optparms{i})
    box off
end

if ishandle(3), close(3); end
figure(3)
plot(Xdata.t, Xdata.F, 'k-', 'linewidth', 1.5); hold on
plot(Xdata.t, IG.Fsei, '-', 'linewidth', 1.5)
plot(out.t, out.F, 'linewidth', 1.5)
plot(t, Fse,':', 'linewidth', 1.5)
box off
xlabel('Time (s)')
ylabel('Force (F_0)')

legend('Target', 'Initial guess', 'Solution', 'Validation', 'location', 'best')


%%
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
    optparms = {'n','kappa', 'kse', 'kse0', 'vmax'};
    
elseif strcmp(model, 'Hill-type no SE')
    optparms = {'n','kappa', 'vmax'};
    
elseif strcmp(model, '4-state XB coop')
    optparms = {'f', 'k11', 'k22', 'k21', 'kon', 'kse', 'kse0', 'kF', 'koop', 'dLcrit', 'Lce0', 'kpe'};
end


end

function[parms] = get_initial_parameters(model, bnds, githubfolder, iF)

[~, modelname] = look_up_model(model);

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Reproduce', 'Process'))
% [~, ~, modelname] = get_model_folder('', mcode, 0);

foldername = fullfile(githubfolder, 'biophysical-muscle-model', 'Reproduce', 'Parameters',fibers{iF});

cd(foldername)
load(['parms_', modelname, '.mat'], 'newparms')

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
    %         parms.J1 = 6.17; parms.koop = 5.7; parms.JF = 1e3; parms.J2 =
    %         200;
    
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
    %         parms.J1 = 6.17; parms.koop = 5.7; parms.JF = 1e3; parms.J2 =
    %         200;
    
    parms.k = 1000;
    parms.b = 3000;
    
    %         parms.k = 500; parms.b = 1000;
    parms.dLcrit = 2;
    
    parms.ps2 = 0;
end
end


%% probably need to make these files (also used in test_model.m)
function[X0] = get_initial_state(model, newparms)

if contains(model, 'XB')
    % assumed initial cross-bridge (XB) state
    Q0 = 1e-3; % fraction of XBs bound
    p0 = 0; % mean strain of bound XBs (power-stroke centered and normalized)
    q0 = .1; % standard deviation strain of bound XBs (power-stroke normalized)

    % find the state in which Fse = Fce + Fpe, given the intial XB state
    [Q00, Q20, lce0, Q10, Fse0, Fpe0, Fce0] = find_steady_state(Q0, p0, q0, newparms, 'regular');
    X0 = [Q00 Q20 Fse0 0 0 0];

else
    X0 = 0;
    
end

if contains(model, '2-state')
    X0(end-1) = 1;
end

end

function[t, Fse, Fpe, Fce] = get_forces_from_state(sol, input, newparms)

    t   = sol.x;
    Fse = sol.y(3,:) * newparms.Fscale;

    dLse = max(newparms.Lse_func(Fse, newparms), 0); % can't be negative
    Lce = interp1(input.t, input.L * newparms.gamma, t) - dLse;

    dLce = Lce - newparms.Lce0;
    Fpe = newparms.kpe * dLce  * newparms.Fscale;
    Fpe(newparms.K*dLce < 10) = newparms.kpe * log(1+exp(dLce(newparms.K*dLce < 10)*newparms.K))/newparms.K  * newparms.Fscale;

    Fce = Fse - Fpe;
end

function[modelfunc, modelname] = look_up_model(model)

if contains(model, 'XB')
    modelfunc = @fiber_dynamics_explicit_length_v2;
else
    modelfunc = @hill;
end

if strcmp(model, '3-state XB coop')
    modelname = 'biophysical_full_regular';

elseif strcmp(model, '2-state XB coop')
    modelname = 'biophysical_thin_regular';

elseif strcmp(model, '2-state XB')
    modelname = 'biophysical_no_regular';

elseif strcmp(model, '4-state XB coop')
    modelname = 'biophysical_full_alternative';
end
end

function[input] = ramp(Ca, T, dt, ISI, A)

    t = 0:dt:T;
    N = length(t);

    dTt = .0383/.4545; % test stretch (= constant)
    dTc = A / .4545; % conditioning stretch

    [tis, Cas, Lis, vis] = create_input(T, dTt, dTc, ISI, Ca, N);

    input.t = tis;
    input.L = Lis;
    input.v = vis;
    input.Cas = Cas;
end
