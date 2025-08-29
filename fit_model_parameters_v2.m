function[parms, out] = fit_model_parameters_v2(opti, optparms, w, data, parms)

% parameters
allparms = {'f','k11','k12','k21','k22','JF','koop','J1','J2', 'kon', 'koff', 'kse','kse0', 'kpe', 'Fpe0'};

% create variables for all parameters
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

% create opti variables for parameters that are fitted
for i = 1:length(optparms)
    eval([optparms{i}, '= opti.variable(1);'])
    eval(['opti.subject_to(',num2str(lb(i)), '<', optparms{i}, '<', num2str(ub(i)),');']);
    eval(['opti.set_initial(',optparms{i},',', num2str(parms.(optparms{i})),');']);
end

% dependent parameter
% J2 = J1 * parms.SRX0 / (1-parms.SRX0);

%% extract input and target
Cas = data.Cas;
vts = data.v;
Fts = data.F;
toc = data.t;
Lts = data.L;
idF = data.idF;
idC = 1:300;

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
Ldi     = interp1(sol.x, sol.y(6,:), toc); % length

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
  
% (slack) controls (defined as above)
dQ0dt  = opti.variable(1,N);
dQ1dt  = opti.variable(1,N); 
dQ2dt  = opti.variable(1,N); 
dNondt = opti.variable(1,N); 
dDRXdt = opti.variable(1,N); 

% extra constraints
opti.subject_to(Q1 - Q0 .* p == 0);
opti.subject_to(Q2 - Q0 .* (p.^2 + q) == 0);

% set initial guess
opti.set_initial(Q0, Q0i);
opti.set_initial(Q1, Q1i);
opti.set_initial(Q2, Q2i);
opti.set_initial(Non, Noni);
opti.set_initial(DRX, DRXi);
opti.set_initial(p, Q1i./Q0i);
opti.set_initial(q, Q2i./Q0i - (Q1i./Q0i).^2);
opti.set_initial(Ld, Ldi);
opti.set_initial(dQ0dt, dQ0dti);
opti.set_initial(dQ1dt, dQ1dti);
opti.set_initial(dQ2dt, dQ2dti);
opti.set_initial(dNondt, dNondti);
opti.set_initial(dDRXdt, dDRXdti);

%% dynamics constraints
F = Q0 + Q1; % cross-bridge force

% specify dynamics using error
error_thin      = ThinEquilibrium(Cas, Q0, Non, dNondt, kon, koff, koop, parms.Noverlap); % thin filament dynamics     
error_thick     = ThickEquilibrium(F, DRX, dDRXdt, J1, J2, JF, parms.Noverlap); % thick filament dynamics
[error_XB, Fdot]= MuscleEquilibrium(Q0, p, q, dQ0dt, dQ1dt, dQ2dt, f, parms.w, k11, k12, k21, k22,  Non, Ld, DRX); % cross-bridge dynamics
error_length    = LengthEquilibrium(Q0, Q1, Fdot, Ld, vts, kse0, kse);

% separate out the XB error
error_Q0 = error_XB(1,:);
error_Q1 = error_XB(2,:);
error_Q2 = error_XB(3,:);

% set errors equal to zero
opti.subject_to(error_thin(:) == 0);
opti.subject_to(error_thick(:) == 0);
opti.subject_to(error_Q0(:) == 0);
opti.subject_to(error_Q1(:) == 0);
opti.subject_to(error_Q2(:) == 0);
opti.subject_to(error_length(:) == 0);

opti.subject_to(J2 - (J1 * parms.SRX0 / (1-parms.SRX0)) == 0);


%% derivative constraints
opti.subject_to((dNondt(1:N-1) + dNondt(2:N))*dt/2 + Non(1:N-1) == Non(2:N));
opti.subject_to((dDRXdt(1:N-1) + dDRXdt(2:N))*dt/2 + DRX(1:N-1) == DRX(2:N));
opti.subject_to((dQ0dt(1:N-1) + dQ0dt(2:N))*dt/2 + Q0(1:N-1) == Q0(2:N));
opti.subject_to((dQ1dt(1:N-1) + dQ1dt(2:N))*dt/2 + Q1(1:N-1) == Q1(2:N));
opti.subject_to((dQ2dt(1:N-1) + dQ2dt(2:N))*dt/2 + Q2(1:N-1) == Q2(2:N));

%% cost
Frel = F * parms.Fscale + kpe * Lts + Fpe0;

% cost function
J = 0;
J = J + w(1) * sum((Frel(idF) - Fts(idF)).^2); % force-velocity fitting
J = J + w(3) * (sum(dQ0dt(idC).^2) + sum(dQ1dt(idC).^2) + sum(dQ2dt(idC).^2)); % regularization term

% optimize
opti.minimize(J); 

%% Solve problem
% options for IPOPT
% options.ipopt.tol = 1*10^(-6);          
% options.ipopt.linear_solver = 'mumps';
% opti.solver('ipopt',options);

% Solve the OCP
p_opts = struct('expand',true);
s_opts = struct('max_iter', 100);
opti.solver('ipopt',p_opts,s_opts);

% visualize 
opti.callback(@(i) plot(toc, [Fts; opti.debug.value(Frel)]))

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
    R.J     = sol.value(J);
    
    % extract the parameters
    for i = 1:length(optparms)
        parms.(optparms{i}) = eval(['sol.value(',optparms{i},');']);
    end

catch
%     sol = opti.debug();

    % Extract the result
    R.Q0    = opti.debug.value(Q0); 
    R.Q1    = opti.debug.value(Q1); 
    R.Q2    = opti.debug.value(Q2); 
    R.dQ0dt = opti.debug.value(dQ0dt); 
    R.dQ1dt = opti.debug.value(dQ1dt); 
    R.dQ2dt = opti.debug.value(dQ2dt); 
    R.F     = opti.debug.value(Frel); 
    R.J     = opti.debug.value(J);
    
    % extract the parameters
    for i = 1:length(optparms)
        parms.(optparms{i}) = eval(['opti.debug.value(',optparms{i},');']);
    end
end

R.Fdot  = R.dQ0dt + R.dQ1dt;
R.t     = 0:dt:(N-1)*dt;

% parms.J2 = sol.value(J2);

%% Test the result: run a forward simulation (sanity check)
[t,x] = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 max(toc)], sol0.y(:,end), xp0, odeopt);
F = (x(:,1) + x(:,2)) * parms.Fscale;
Fn = interp1(t, F, toc) + parms.kpe * Lts + parms.Fpe0;

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
