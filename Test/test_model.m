clear all; close all; clc
cd(fileparts(which('test_model.m')));
cd ..
addpath(genpath('Common'))
repofolder = cd;

%% step 1: specify inputs
% the biophysical models have 3 inputs, namely:
% 1) calcium concentration
% 2) fiber length
% 3) fiber velocity
% note: the latter two depend on each other, because
% fiber velocity is the time-derivative of fiber length
%
% in the following example, three phases of input trajectories are used:
% 1) typical protocol from our manuscript
% 2) an isometric phase in which the fiber is held at constant length
% 3) a sinusoidal phase in which the fiber undergoes sinusoidal length
% changes

% applicable to all phases
pCa     = 9;          % assumed constant and applies to both phases
Ca      = 10^(6-pCa);   % (uM)
dt      = 1/1000;       % sample time (s)
T       = 1;            % duration (s)
L0      = 0;

% phase-specific
A   = .2;              % length change amplitude (L0)
RT  = .1;             % recovery time (s)
f   = 2;              % sinusoidal frequency (Hz)
v   = 1;            % velocity (L0/s)

% specify all input phases
input(1) = ramp(Ca, T, dt, RT, A, v);
input(2) = isometric(Ca, T, dt, L0);
input(3) = sinusoidal(Ca, T, dt, f, A, L0);

% show length and velocity traces
t0 = 0;
color = lines(3);
for i = 1:length(input) % loop over phases
    if i > 1, t0 = input(i-1).t(end) + t0;
    end
    
    figure(1)
    subplot(411)
    plot(input(i).t+t0, input(i).Cas, 'color', color(1,:), 'linewidth', 1.5); hold on; box off
    ylabel('[Ca]^{2+} (\muM)')
    xline(t0,'k--')
    title('Calcium concentration')
    
    subplot(412)
    plot(input(i).t+t0, input(i).L, 'color', color(1,:), 'linewidth', 1.5); hold on; box off
    ylabel('\DeltaFiber length (L_0)')
    xline(t0,'k--')
    title('Fiber length change')
    
    subplot(413)
    plot(input(i).t+t0, input(i).v, 'color', color(1,:), 'linewidth', 1.5); hold on; box off
    ylabel('Fiber velocity (L_0/s)')
    xline(t0,'k--')
    title('Fiber velocity')
    
end

% note: this is just an example. you are free to try any input that you
% like, such as:
% - different calcium concentration
% - different amplitude and frequency of sinusoidal length changes
%
% but also:
% - variable calcium concentration
% - non-sinusoidal fiber length changes
% 
% if you use a custom input, make sure that you have a struct (array) called
% 'input' that contains the following fields:
%     t: [1×N double]
%     L: [1×N double]
%     v: [1×N double]
%    Ca: [1×N double]
% here, N is the number of data points, t should be increasing
% equidistantly (i.e. diff(t) = constant), and v should be the time
% derivtive of L with respect to t.

%% step 2: specify model function
% you can choose between the following models:
% Hill-type SE, Hill-type no SE, 2-state XB, 2-state XB coop, 3-state XB coop, 4-state XB coop
model   = '3-state XB coop'; % see options above
odetype = 'explicit'; % type of differential equations
method  = 'approximated'; % solution method
% method = 'discretized'; % solution method
[modelfunc, odefunc, modelname] = look_up_model(model, odetype, method);

%% step 3: specify model parameters
% here we look up the parameters for a given fiber and model
cd(fullfile(repofolder, 'Reproduce', 'Parameters'))
% fibername       = uigetdir('Pick a fiber');
fibername = '7Aug2018a';
fullfilename    = fullfile(fibername, ['parms_', modelname, '.mat']);
load(fullfilename, 'redparms')
newparms = complete_parms(redparms);

%% step 4: simulate model
[sol, out] = simulate_model(model, odefunc, modelfunc, input, newparms);

% visualize
figure(1)
subplot(414)
t0 = 0;
for i = 1:length(input)
    if i > 1, t0 = input(i-1).t(end) + t0;
    end

    % plot the forces
    plot(out(i).t+t0, out(i).F, 'color', color(1,:), 'linewidth', 1.5); hold on; box off
end

xlabel('Time (s)')
ylabel('Force (F_0)')
xline(t0,'k--')
title('Fiber force')

