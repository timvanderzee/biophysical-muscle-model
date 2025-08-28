function[parms, out] = fit_model_parameters_v2(opti, optparms, w, data, parms)

% define desired force-velocity properties
% vmax = 5; % maximal shortening velocity [L0/s]

% parms.J2 = 200;

%% parameters
% parms.Fscale = 2;

% parameters
allparms = {'f','k11','k12','k21','k22','JF','koop','J1','J2', 'kon', 'koff', 'kse','kse0', 'kpe', 'Fpe0'};

for i = 1:length(allparms)
    eval([allparms{i}, ' = ', num2str(parms.(allparms{i})),';'])
end

% bounds on parameter values
lb = .1 * ones(1, length(optparms));
ub = 2e3 * ones(1, length(optparms));

% exception for exponential coefficient
for i = 1:length(optparms)
    if strcmp(optparms{i}, 'k22') || strcmp(optparms{i}, 'k12') || strcmp(optparms{i}, 'kse') || strcmp(optparms{i}, 'kse0')
        lb(i) = 0;
        ub(i) = 5;
    end
end

for i = 1:length(optparms)
    if strcmp(optparms{i}, 'kpe') || strcmp(optparms{i}, 'Fpe0')
        lb(i) = -1;
        ub(i) = 5;
    end
end

for i = 1:length(optparms)
    eval([optparms{i}, '= opti.variable(1);'])
    eval(['opti.subject_to(',num2str(lb(i)), '<', optparms{i}, '<', num2str(ub(i)),');']);
    eval(['opti.set_initial(',optparms{i},',', num2str(parms.(optparms{i})),');']);
end

% enforce constraint on J2
% SRX0 = .8;
% J2 = J1 * SRX0 / (1-SRX0);

%% design velocity input vector
% this is for testing both force-velocity and history-dependent properties
% N = 1000; % number of nodes 
% [vts, Fts, toc, idF, idFd] = design_length_input_vector(vmax, RT, V_rel, N);
% dt = mean(diff(toc));

Cas = data.Cas;
vts = data.v;
Fts = data.F;
toc = data.t;
Lts = data.L;

idF = data.idF;
idC = data.idC;

N = length(toc);
dt = mean(diff(toc));

%% obtain initial guess
% intial guess is obtained through running a forward simulation with the
% first, simulate an isometric contraction
parms.vts = [0 0];
parms.ti = [0 1];
parms.Cas = Cas(1) * [1 1];

x0 = 1e-3 * ones(6,1);
xp0 = zeros(size(x0));
odeopt = odeset('maxstep', 3e-3);
sol0 = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 1], x0, xp0, odeopt);
% [~,xdot0] = deval(sol0, sol0.x);

% next, simulate response to specified velocity input vector
parms.vts = vts;
parms.ti = toc;
parms.Cas = Cas;

sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 max(toc)], sol0.y(:,end), xp0, odeopt);
[~,xdot] = deval(sol, sol.x);

% interpolate solution to time nodes
Q0i     = interp1(sol.x, sol.y(1,:), toc); % zero-order moment
Q1i     = interp1(sol.x, sol.y(2,:), toc); % first-order moment
Q2i     = interp1(sol.x, sol.y(3,:), toc); % second-order moment
Noni    = interp1(sol.x, sol.y(4,:), toc); % thin filament activation
DRXi    = interp1(sol.x, sol.y(5,:), toc); % thick filament activation
Ldi    = interp1(sol.x, sol.y(6,:), toc); % length

dQ0dti  = interp1(sol.x, xdot(1,:), toc); % zero-order moment time derivative
dQ1dti  = interp1(sol.x, xdot(2,:), toc); % first-order moment time derivative
dQ2dti  = interp1(sol.x, xdot(3,:), toc); % second-order moment time derivative
dNondti = interp1(sol.x, xdot(4,:), toc); % thin filament activation time derivative
dDRXdti = interp1(sol.x, xdot(5,:), toc); % thick filament activation time derivative

Fi = (Q0i + Q1i) * parms.Fscale;
Fidot = (dQ0dti + dQ1dti) * parms.Fscale;

%% Fit cross-bridge rates using direct collocation
% define opti states (defined as above)
Q0  = opti.variable(1,N);
Q1  = opti.variable(1,N); 
Q2  = opti.variable(1,N); 
Non = opti.variable(1,N);
DRX = opti.variable(1,N);
Ld  = opti.variable(1,N);

% define extra variables
p  = opti.variable(1,N); % mean strain of the distribution
q  = opti.variable(1,N); % standard deviation strain of the distribution
% F  = opti.variable(1,N); % force of the distribution
  
% (Slack) controls (defined as above)
dQ0dt  = opti.variable(1,N);
dQ1dt  = opti.variable(1,N); 
dQ2dt  = opti.variable(1,N); 
dNondt = opti.variable(1,N); 
dDRXdt = opti.variable(1,N); 

% Inequality constraints
% opti.subject_to(Q0 >= 0);
% opti.subject_to(Q1 >= -Q0);
% opti.subject_to(q >= 0);
% opti.subject_to(Non >= 0);
% opti.subject_to(Non <= 1);
% opti.subject_to(DRX >= 0);
% opti.subject_to(DRX <= 1);

% start from rest
% opti.subject_to(Q0(1) == x0(1));
% opti.subject_to(Q1(1) == x0(2));
% opti.subject_to(Q2(1) == x0(3));
% opti.subject_to(Non(1) == x0(4));
% opti.subject_to(DRX(1) == x0(5));

% Extra constraints
% eps = 0;
opti.subject_to(Q1 - Q0 .* p == 0);
opti.subject_to(Q2 - Q0 .* (p.^2 + q) == 0);
% opti.subject_to(F - parms.Fscale * (Q0 + Q1) == 0);

% Set initial guess states based on forward simulation results
% Q0i = 0.5 * ones(1,N);
% Q1i = 0.5 * ones(1,N);
% Q2i = 0.6 * ones(1,N);
% Noni = 0.6 * ones(1,N);
% DRXi = 0.6 * ones(1,N);
% 
% dQ0dti = 0 * ones(1,N);
% dQ1dti = 0 * ones(1,N);
% dQ2dti = 0 * ones(1,N);
% dNondti = 0 * ones(1,N);
% dDRXdti = 0 * ones(1,N);

opti.set_initial(Q0, Q0i);
opti.set_initial(Q1, Q1i);
opti.set_initial(Q2, Q2i);
opti.set_initial(Non, Noni);
opti.set_initial(DRX, DRXi);
opti.set_initial(p, Q1i./Q0i);
opti.set_initial(q, Q2i./Q0i - (Q1i./Q0i).^2);
opti.set_initial(Ld, Ldi);
% opti.set_initial(F, parms.Fscale * (Q0i + Q1i));

opti.set_initial(dQ0dt, dQ0dti);
opti.set_initial(dQ1dt, dQ1dti);
opti.set_initial(dQ2dt, dQ2dti);
opti.set_initial(dNondt, dNondti);
opti.set_initial(dDRXdt, dDRXdti);

%% dynamics constraints
F = Q0 + Q1; % cross-bridge force

% k = 100;
% Nonc = log(1+exp(Non*k))/k;
% DRXc = log(1+exp(DRX*k))/k;
% Q0c = log(1+exp(Q0*k))/k;
% Fc = log(1+exp(F*k))/k;
% qc = log(1+exp(q*k))/k;

Nonc = Non;
DRXc = DRX;
Q0c = Q0;
Fc = F;
qc = q;

% error = [];
error_thin      = ThinEquilibrium(Cas, Q0c, Nonc, dNondt, kon, koff, koop, parms.Noverlap); % thin filament dynamics     
error_thick     = ThickEquilibrium(Fc, DRXc, dDRXdt, J1, J2, JF, parms.Noverlap); % thick filament dynamics
[error1, Fdot]  = MuscleEquilibrium(Q0c, p, qc, dQ0dt, dQ1dt, dQ2dt, f, parms.w, k11, k12, k21, k22,  Nonc, Ld, DRXc); % cross-bridge dynamics
error_length    = LengthEquilibrium(Q0c, Q1, Fdot, Ld, vts, kse0, kse);

error_Q0 = error1(1,:);
error_Q1 = error1(2,:);
error_Q2 = error1(3,:);
% error       = [error; error_thin(:); error_thick(:); error1(:)];
% opti.subject_to(error == 0);

opti.subject_to(error_thin(:) == 0);
opti.subject_to(error_thick(:) == 0);
opti.subject_to(error_Q0(:) == 0);
opti.subject_to(error_Q1(:) == 0);
opti.subject_to(error_Q2(:) == 0);
opti.subject_to(error_length(:) == 0);

%% derivative constraints
opti.subject_to((dNondt(1:N-1) + dNondt(2:N))*dt/2 + Non(1:N-1) == Non(2:N));
opti.subject_to((dDRXdt(1:N-1) + dDRXdt(2:N))*dt/2 + DRX(1:N-1) == DRX(2:N));
opti.subject_to((dQ0dt(1:N-1) + dQ0dt(2:N))*dt/2 + Q0(1:N-1) == Q0(2:N));
opti.subject_to((dQ1dt(1:N-1) + dQ1dt(2:N))*dt/2 + Q1(1:N-1) == Q1(2:N));
opti.subject_to((dQ2dt(1:N-1) + dQ2dt(2:N))*dt/2 + Q2(1:N-1) == Q2(2:N));

%% cost
Frel = Fc * parms.Fscale + kpe * Lts + Fpe0;
% Freldot = dQ0dt + dQ1dt;

% idF = 1:length(Frel);

% idF = find((toc > 0.2 & toc < 0.8) | (toc > 1.2 & toc < 1.8) | (toc > 2.2 & toc < 2.8)); 
% idF = find((toc > .5 & toc < 1.5) | (toc > 3.5 & toc < 4.5) | (toc > 6.5 & toc < 7.5)); 

idC = 1:300;

% cost function
J = 0;
J = J + w(1) * sum((Frel(idF) - Fts(idF)).^2); % force-velocity fitting
% J = J + w(2) * sum((SRS_rel - Freldot(idFd(2))/Freldot(idFd(1))).^2); % history-dependence fitting
J = J + w(3) * (sum(dQ0dt(idC).^2) + sum(dQ1dt(idC).^2) + sum(dQ2dt(idC).^2)); % regularization term
% J = J + .01 * f.^2;

opti.minimize(J); 

%% Solve problem
% options for IPOPT
options.ipopt.tol = 1*10^(-6);          
options.ipopt.linear_solver = 'mumps';
% opti.solver('ipopt',options);

% Solve the OCP
p_opts = struct('expand',true);
s_opts = struct('max_iter', 100);
opti.solver('ipopt',p_opts,s_opts);

% visualize 
opti.callback(@(i) plot(toc, [Fts; opti.debug.value(Frel)]))

% opti.callback(@(i) plot((opti.debug.value(qc)))); hold on

try
    sol = opti.solve();  

    % Extract the result
    R.Q0    = sol.value(Q0); 
    R.Q1    = sol.value(Q1); 
    R.Q2    = sol.value(Q2); 
    R.dQ0dt = sol.value(dQ0dt); 
    R.dQ1dt = sol.value(dQ1dt); 
    R.dQ2dt = sol.value(dQ2dt); 
    R.F     = sol.value(Frel); 
    R.Fdot  = R.dQ0dt + R.dQ1dt;
    R.t     = 0:dt:(N-1)*dt;
    
    R.J = sol.value(J);

catch
    sol = opti.debug();

    % Extract the result
    R.Q0    = opti.debug.value(Q0); 
    R.Q1    = opti.debug.value(Q1); 
    R.Q2    = opti.debug.value(Q2); 
    R.dQ0dt = opti.debug.value(dQ0dt); 
    R.dQ1dt = opti.debug.value(dQ1dt); 
    R.dQ2dt = opti.debug.value(dQ2dt); 
    R.F     = opti.debug.value(Frel); 
    R.Fdot  = R.dQ0dt + R.dQ1dt;
    R.t     = 0:dt:(N-1)*dt;

    R.J = opti.debug.value(J);
end

R.Non = sol.value(Non);
R.dNondt = sol.value(dNondt);
R.error_thin = sol.value(error_thin);

%% Test the result
% extract the parameters
for i = 1:length(optparms)
    parms.(optparms{i}) = eval(['sol.value(',optparms{i},');']);
end


%% Test errors
% k = 20;
% Nonc = log(1+exp(R.Non*k))/k;
% DRXc = log(1+exp(R.DRX*k))/k;
% Q0c = log(1+exp(R.Q0*k))/k;
% Fc = log(1+exp(R.F*k))/k;

% error_thin = ThinEquilibrium(Cas, R.Q0, R.Non, R.dNondt, parms.kon, parms.koff, koop, parms.Noverlap);
% error_thin2 = ThinEquilibrium(Cas, Q0c, Nonc, R.dNondt, parms.kon, parms.koff, koop, parms.Noverlap);

% run a forward simulation
osol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 max(toc)], sol0.y(:,end), xp0, odeopt);
t = osol.x;
x = osol.y;
F = (x(1,:) + x(2,:)) * parms.Fscale;
Fn = interp1(t, F, toc) + parms.kpe * Lts + parms.Fpe0;

% [~,xdot] = deval(osol, t);
% Fdot = xdot(1,:) + xdot(2,:);

%% Visualize
subplot(311)
plot(R.t, vts, 'linewidth',1.5); hold on
ylabel('Velocity')
title('Velocity')

subplot(312)
plot(R.t(idF), Fts(idF), 'k.', 'markersize',10); hold on 
plot(toc, Fi, ':', 'linewidth',1.5); 
plot(R.t, R.F, 'linewidth',1.5); hold on
plot(toc, Fn, ':', 'linewidth',1.5); 
ylabel('Force')
title('Force')
legend('Target','Initial guess','Result','Simulated','location','best')
legend boxoff

subplot(313);
plot(toc, Fidot); hold on
plot(R.t, R.Fdot*parms.Fscale, 'linewidth',1.5);
ylabel('Force-rate')
title('Force-rate')

for i = 1:3
    subplot(3,1,i); 
    hold on
    box off
    xlabel('Time (s)')
    xlim([0 max(toc)])
end

set(gcf,'units','normalized','position', [.1 .1 .4 .8])

%% output
out.v = vts(idF);
out.F = R.F(idF);
out.Ft = Fts(idF);
out.cost = R.J;

end
