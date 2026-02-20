clear all; close all; clc
cd ..
addpath(genpath('Common'))
cd ..
githubfolder = cd;

%% step 1: specify which model we want to fit
model = '3-state XB coop'; 
[modelfunc, modelname] = look_up_model(model);

%% step 2: specify which parameters we are fitting, and which bounds we are using
optparms    = {'f', 'k11', 'J1', 'kse'};
bnds.f      = [1e0 1e3];
bnds.k11    = [1e0 1e3];
bnds.J1     = [1e0 1e3];
bnds.kse    = [1e-3 1e1];

%% step 3: obtain (initial) parameter values
iF = 2; % fiber number (allowed: [2,3,5,6,7,8,11]);
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
foldername = fullfile(githubfolder, 'biophysical-muscle-model', 'Reproduce', 'Parameters',fibers{iF});

cd(foldername)
load(['parms_', modelname, '.mat'], 'newparms')
oldparms = newparms;

%% step 4: specify which data we want to fit on
% similar to as in test_model.m, we need an input struct called
% 'input' that contains the following fields:
%     t: [1×N double]
%     L: [1×N double]
%     v: [1×N double]
%    Ca: [1×N double]
%
% in addition to these fields, we need the following field:
%    F: [1×N double]
%
% this last field contains the force that we want to fit. 

% as in test_model.m, we can create inputs
T = 1;          % total duration (s)
dt = 1/1000;    % sample time (s)
Ca = 1;         % calcium concentration (uM)
RT = .1;        % recovery time(s)
A = .04;        % stretch amplitude (L0)

% make an input struct
input = ramp(Ca, T, dt, RT, A);
color = lines(3);

% at this time, we still need to add the force data that we'd like to fit
% here, we will fit an adjusted simulated force trajectory. we start with
% the simulated force corresponding to the above input (input), the model
% (model), and the model parameters (oldparms). next, we adjust this force
% through multiplying it with a scaling factor, and adding an offset and
% noise. as a result of these adjustments, there will be a mismatch between
% the initial model force and the target model force. this mismatch can be
% reduced by re-fitting the model parameters
[sol, out] = simulate_model(model, modelfunc, input, oldparms);

% adjustment parameters
offset = -.2;
scaling = 1.2;
A_noise = .01;

% create adjust target force
input.F = offset + interp1(sol.x, out.F, input.t) * scaling + randn(size(input.t)) * A_noise;

% visualize data to be fitted
if ishandle(1), close(1); end
figure(1)
subplot(411)
plot(input.t, input.Cas, 'color', color(1,:), 'linewidth', 1.5); hold on; box off
ylabel('[Ca]^{2+} (\muM)')
title('Calcium concentration')

subplot(412)
plot(input.t, input.L, 'color', color(1,:), 'linewidth', 1.5); hold on; box off
ylabel('\DeltaFiber length (L_0)')
title('Fiber length change')

subplot(413)
plot(input.t, input.v, 'color', color(1,:), 'linewidth', 1.5); hold on; box off
ylabel('Fiber velocity (L_0/s)')
title('Fiber velocity')

subplot(414)
plot(input.t, input.F, 'color', color(1,:), 'linewidth', 1.5); hold on; box off
ylabel('Fiber force (F_0)')
title('Fiber force')
    
%% step 5: obtain initial guess for model states
cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Fitting'))
IG = get_initial_guess(model, modelfunc, input, oldparms, 0);

subplot(414)
plot(input.t, IG.Fsei*oldparms.Fscale, 'color', color(2,:), 'linewidth', 1.5); hold on; box off
legend('Target', 'Initial guess', 'location', 'best')

%% step 6: do fitting
% define weight vector
weights = 100;

% fit model parameters
[newparms, out, opti] = fit_model_parameters_v2(model, oldparms, optparms, bnds, input, IG, weights);

%% step 7: validate through simulating with the obtained parameters
[sol, sim_out] = simulate_model(model, modelfunc, input, newparms);

%% visualize output
if ishandle(2), close(2); end
figure(2)
for i = 1:length(optparms)
    nexttile
    bar(categorical({'old', 'new'}), [oldparms.(optparms{i}) newparms.(optparms{i})])
    title(optparms{i})
    box off
end

if ishandle(3), close(3); end
figure(3)
plot(input.t, input.F, 'k-', 'linewidth', 1.5); hold on
plot(input.t, IG.Fsei * oldparms.Fscale, '-', 'linewidth', 1.5)
plot(out.t, out.F, 'linewidth', 1.5)
plot(sol.x, sim_out.F,':', 'linewidth', 1.5)
box off
xlabel('Time (s)')
ylabel('Force (F_0)')

legend('Target', 'Initial guess', 'Solution', 'Validation', 'location', 'best')
