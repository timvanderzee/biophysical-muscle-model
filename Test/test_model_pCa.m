clear all; close all; clc
cd(fileparts(which('test_model.m')));
cd ..
addpath(genpath('Common'))
repofolder = cd;

%% step 1: specify inputs
% applicable to all phases
pCas     = flip([9 7:-.1:5 4.5]);          % assumed constant and applies to both phases
Cas      = 10.^(6-pCas);   % (uM)
dt      = 1/1000;       % sample time (s)
T       = 1;            % duration (s)
L0      = 0;

% phase-specific
A   = .2;              % length change amplitude (L0)
RT  = .1;             % recovery time (s)
f   = 2;              % sinusoidal frequency (Hz)
v   = 1;            % velocity (L0/s)

% specify all input phases
for i = 1:length(Cas)
    input(i) = isometric(Cas(i), T, dt, L0);
end

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

%% step 2: specify model function
% you can choose between the following models:
% Hill-type SE, Hill-type no SE, 2-state XB, 2-state XB coop, 3-state XB coop, 4-state XB coop
model   = '3-state XB coop'; % see options above
odetype = 'explicit'; % type of differential equations
method  = 'approximated'; % solution method
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
tic
[sol, out] = simulate_model(model, odefunc, modelfunc, input, newparms);
toc

% visualize
figure(1)
subplot(414)
t0 = 0;
for i = 1:length(input)
    if i > 1, t0 = input(i-1).t(end) + t0;
    end

    % plot the forces
    plot(out(i).t+t0, out(i).F, 'color', color(1,:), 'linewidth', 1.5); hold on; box off
    xline(t0,'k--')

    Fss(i) = out(i).Fce(end);
end

xlabel('Time (s)')
ylabel('Force (F_0)')
xline(t0,'k--')
title('Fiber force')

figure(2)
plot(pCas, Fss, 'o-', 'markerfacecolor', [1 1 1])
xlabel('pCa')
ylabel('Isometric force')
box off

