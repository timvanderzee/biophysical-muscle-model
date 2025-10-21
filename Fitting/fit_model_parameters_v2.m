function[parms, out] = fit_model_parameters_v2(opti, optparms, w, data, parms, IG)

% parameters
allparms = {'f','k11','k12','k21','k22','JF','koop','J1','J2', 'kon', 'koff', 'kse','kse0', 'kpe', 'Fpe0','b','k','dLcrit', 'gamma'};

% create variables for all parameters
for i = 1:length(allparms)
    eval([allparms{i}, ' = ', num2str(parms.(allparms{i})),';'])
end

% % bounds on parameter values
% lb = .1 * ones(1, length(optparms));
% ub = 1e3 * ones(1, length(optparms));
% 
% % exception for exponential coefficient
% for i = 1:length(optparms)
%     if strcmp(optparms{i}, 'k22') || strcmp(optparms{i}, 'k12') || strcmp(optparms{i}, 'kse') || strcmp(optparms{i}, 'kse0')
%         lb(i) = 1e-4;
%         ub(i) = 5;
%     end
% end

for i = 1:length(optparms)
    lb(i) = .2 * parms.(optparms{i});
    ub(i) = 5 * parms.(optparms{i});
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
idC = [1:300 data.idC];

N = length(toc);
dt = mean(diff(toc));

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
F  = opti.variable(1,N);
F0dot  = opti.variable(1,N);

% (slack) controls (defined as above)
dQ0dt  = opti.variable(1,N);
dQ1dt  = opti.variable(1,N); 
dQ2dt  = opti.variable(1,N); 
dNondt = opti.variable(1,N); 
dDRXdt = opti.variable(1,N); 

% extra constraints
opti.subject_to(dQ0dt + dQ1dt - F0dot - Ld .* Q0 == 0);
opti.subject_to(Q0 + Q1 - F == 0);
opti.subject_to(Q1 - Q0 .* p == 0);
opti.subject_to(Q2 - Q0 .* (p.^2 + q) == 0);

% opti.subject_to(p >= -5);
% opti.subject_to(p <= 5)

opti.subject_to(q >= 0);
opti.subject_to(Q0 >= 0);
opti.subject_to(F >= 0);
opti.subject_to(Non >= 0);
opti.subject_to(DRX >= 0);

% opti.subject_to(Q0 <= 1);
% opti.subject_to(Non <= 1);
% opti.subject_to(DRX <= 1);


% set initial guess
opti.set_initial(F, IG.Fi);
opti.set_initial(Q0, IG.Q00i);
opti.set_initial(Q1, IG.Q1i);
opti.set_initial(Q2, IG.Q2i);
opti.set_initial(Non, IG.Noni);
opti.set_initial(DRX, IG.DRXi);
opti.set_initial(p, IG.pi);
opti.set_initial(q, IG.qi);
opti.set_initial(Ld, IG.Ldi);
opti.set_initial(dQ0dt, IG.dQ0dti);
opti.set_initial(dQ1dt, IG.dQ1dti);
opti.set_initial(dQ2dt, IG.dQ2dti);
opti.set_initial(dNondt, IG.dNondti);
opti.set_initial(dDRXdt, IG.dDRXdti);

opti.set_initial(F0dot, IG.F0doti);

if parms.b > 0
    R = opti.variable(1,N);
    dRdt = opti.variable(1,N);

%     opti.subject_to(R >= 0);
    opti.set_initial(R, IG.Ri);
    opti.set_initial(dRdt, IG.dRdti);
else
    R = zeros(1, N);
    dRdt = zeros(1,N);
end

%% dynamics constraints
% F = Q0 + Q1; % cross-bridge force

% specify dynamics using error
error_thin      = ThinEquilibrium(Cas, Q0, Non, dNondt, kon, koff, koop, parms.Noverlap); % thin filament dynamics     
error_thick     = ThickEquilibrium(Q0, dQ0dt, F, DRX, dDRXdt, J1, J2, JF, parms.Noverlap); % thick filament dynamics
[error_Q0, error_Q1, error_Q2, error_R] = MuscleEquilibrium(Q0, Q1, p, q, dQ0dt, dQ1dt, dQ2dt, f, parms.w, k11, k12, k21, k22,  Non, Ld, DRX, dRdt, b, k, R, dLcrit); % cross-bridge dynamics
error_length    = LengthEquilibrium(Q0, F, F0dot, Ld, vts, kse0, kse, parms.gamma);

% set errors equal to zero
opti.subject_to(error_thin(:) == 0);
opti.subject_to(error_thick(:) == 0);
opti.subject_to(error_Q0(:) == 0);
opti.subject_to(error_Q1(:) == 0);
opti.subject_to(error_Q2(:) == 0);
opti.subject_to(error_length(:) == 0);

if parms.b > 0
   opti.subject_to(error_R(:) == 0);
end

% opti.subject_to(J2 - J1 * parms.SRX0 / (1-parms.SRX0) == 0);

%% derivative constraints
opti.subject_to((dNondt(1:N-1) + dNondt(2:N))*dt/2 + Non(1:N-1) == Non(2:N));
opti.subject_to((dDRXdt(1:N-1) + dDRXdt(2:N))*dt/2 + DRX(1:N-1) == DRX(2:N));
opti.subject_to((dQ0dt(1:N-1) + dQ0dt(2:N))*dt/2 + Q0(1:N-1) == Q0(2:N));
opti.subject_to((dQ1dt(1:N-1) + dQ1dt(2:N))*dt/2 + Q1(1:N-1) == Q1(2:N));
opti.subject_to((dQ2dt(1:N-1) + dQ2dt(2:N))*dt/2 + Q2(1:N-1) == Q2(2:N));

if parms.b > 0
    opti.subject_to((dRdt(1:N-1) + dRdt(2:N))*dt/2 + R(1:N-1) == R(2:N));
end

%% cost
Frel = F * parms.Fscale + kpe * Lts + Fpe0;
Fcost = (Frel(idF) - Fts(idF)).^2;

% cost function
J = 0;
J = J + w(1) * sum(Fcost); % force-velocity fitting
J = J + w(3) * (sum(dQ0dt(idC).^2) + sum(dQ1dt(idC).^2) + sum(dQ2dt(idC).^2)); % regularization term

% optimize
opti.minimize(J); 

%% Solve problem
% options for IPOPT
% options.ipopt.tol = 1*10^(-6);          
% options.ipopt.linear_solver = 'mumps';
% opti.solver('ipopt',options);

% Solve the OCP
p_opts = struct('detect_simple_bounds', true);
s_opts = struct('max_iter', 100);
opti.solver('ipopt',p_opts,s_opts);

% visualize 
opti.callback(@(i) plot(toc, [Fts; opti.debug.value(Frel)]))

try
    sol = opti.solve();  

    % Extract the result
    out.Q0    = sol.value(Q0); 
    out.Q1    = sol.value(Q1); 
    out.Q2    = sol.value(Q2); 
    out.dQ0dt = sol.value(dQ0dt); 
    out.dQ1dt = sol.value(dQ1dt); 
    out.dQ2dt = sol.value(dQ2dt); 
    out.F     = sol.value(Frel); 
    out.J     = sol.value(J);
    out.Fcost = sol.value(Fcost);
    
    % extract the parameters
    for i = 1:length(optparms)
        parms.(optparms{i}) = eval(['sol.value(',optparms{i},');']);
    end

catch
%     sol = opti.debug();

    % Extract the result
    out.Q0    = opti.debug.value(Q0); 
    out.Q1    = opti.debug.value(Q1); 
    out.Q2    = opti.debug.value(Q2); 
    out.dQ0dt = opti.debug.value(dQ0dt); 
    out.dQ1dt = opti.debug.value(dQ1dt); 
    out.dQ2dt = opti.debug.value(dQ2dt); 
    out.F     = opti.debug.value(Frel); 
    out.J     = opti.debug.value(J);
    out.Fcost = opti.debug.value(Fcost);
    
    % extract the parameters
    for i = 1:length(optparms)
        parms.(optparms{i}) = eval(['opti.debug.value(',optparms{i},');']);
    end
end

out.Fdot  = out.dQ0dt + out.dQ1dt;
out.t     = 0:dt:(N-1)*dt;

% parms.J2 = sol.value(J2);


%% Visualize
subplot(311)
plot(out.t, vts, 'linewidth',1.5); hold on
ylabel('Velocity')
title('Velocity')

subplot(312)
plot(out.t(idF), Fts(idF), 'k.', 'markersize',10); hold on 
% plot(toc, Fi * parms.Fscale + parms.kpe * Lts + parms.Fpe0, ':', 'linewidth',1.5); 
plot(out.t, out.F, 'linewidth',1.5); hold on
% plot(toc, Fn, ':', 'linewidth',1.5); 
ylabel('Force')
title('Force')
% legend('Target','Initial guess','Result','Simulated','location','best')
% legend boxoff

subplot(313);
% plot(toc, Fdoti*parms.Fscale); hold on
plot(out.t, out.Fdot*parms.Fscale, 'linewidth',1.5);
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
% out.v = vts(idF);
% out.F = R.F(idF);
% out.Ft = Fts(idF);
% out.cost = R.J;

end
