clear all; close all; clc
cd ..
addpath(genpath('Functions'))
addpath(genpath('Common'))

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
pCa = 6.1;          % assumed constant and applies to both phases
Ca = 10^(6-pCa);    % (uM)
dt = 1/1000;        % sample time (s)
T = 2;              % duration (s)
L0 = 0;

% phase-specific
A = .04;            % length change amplitude (L0)
RT = .1;            % recovery time (s)
f = 2;              % sinusoidal frequency (Hz)

% specify all input phases
input(1) = ramp(Ca, T, dt, RT, A);
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
    plot(input(i).t+t0, input(i).Ca, 'color', color(1,:), 'linewidth', 1.5); hold on; box off
    ylabel('[Ca]^{2+} (\muM)')
    xline(t0,'k--')
    title('Calcium concentration')
    
    subplot(412)
    plot(input(i).t+t0, input(i).L, 'color', color(1,:), 'linewidth', 1.5); hold on; box off
    ylabel('Fiber velocity (L_0)')
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
% 2-state XB, 2-state XB coop, 3-state XB coop, 4-state XB coop
model = '3-state XB coop'; % see options above
[modelfunc, modelname] = look_up_model(model);

%% step 3: specify model parameters
% here we look up the parameters for a given fiber and model
fiber = '7Aug2018a';
githubfolder = cd;
load(fullfile(cd, 'Parameters', fiber, ['parms_', modelname, '.mat']))

%% step 4: determine initial state
% we need to define the model states at t = 0
X0 = get_initial_state(model, newparms);
sol = ode15s(@(t,y,yp) modelfunc(t,y, newparms), [0 .001], X0, []); % just to define sol

%% step 5: simulate model
for i = 1:length(input) % loop over phases
    
    % obtain the inputs for this phase
    newparms.ti     = input(i).t;
    newparms.Lts    = input(i).L * newparms.gamma;
    newparms.vts    = input(i).v;
    newparms.Cas    = input(i).Ca;
    
    % simulate
    sol(i) = ode15s(@(t,y,yp) modelfunc(t,y, newparms), [0 max(newparms.ti)], X0, []);
    
    % save final state (used as initial state for next simulation)
    X0 = sol(i).y(:,end);
end

%% step 6: calculate model forces from simulation output
t0 = 0;
for i = 1:length(input)
    if i > 1, t0 = input(i-1).t(end) + t0;
    end
    
    % obtain forces
    [t, Fse, Fpe, Fce] = get_forces_from_state(sol(i), input(i), newparms);
    
    % plot the forces
    figure(1)
    subplot(414)
    plot(t+t0, Fse, 'color', color(1,:), 'linewidth', 1.5); hold on; box off
%     plot(t+t0, Fpe, 'color', color(2,:), 'linewidth', 1.5);
%     plot(t+t0, Fce, 'color', color(3,:), 'linewidth', 1.5);
    
    xlabel('Time (s)')
    ylabel('Force (F_0)')
    xline(t0,'k--')
    title('Fiber force')
end

% legend('F_{SE}', 'F_{PE}', 'F_{CE}', 'location', 'best')

%% functions
function[t, Fse, Fpe, Fce] = get_forces_from_state(sol, input, newparms)

    t   = sol.x;
    Fse = sol.y(3,:) * newparms.Fscale;

    dLse = max(newparms.Lse_func(Fse, newparms), 0); % can't be negative
    Lce = interp1(input.t, input.L, t) - dLse;

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

function[input] = isometric(Ca, T, dt, L0)

    input.t  = 0:dt:T;
    input.L  = zeros(size(input.t))  + L0;
    input.v  = zeros(size(input.t));
    input.Ca = ones(size(input.t)) * Ca;

end

function[input] = sinusoidal(Ca, T, dt, f, A, L0)

    input.t = 0:dt:T;
    input.L = A * (-.5 * cos(2*pi*f*input.t) + .5) + L0;
    input.v = A * 2*pi*f*.5*sin(2*pi*f*input.t);
    input.Ca = ones(size(input.t)) * Ca;

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
    input.Ca = Cas;
end

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

