function[parms, out] = fit_model_parameters_v2(optparms, w, data, parms, IG, bnds)

import casadi.*

% initialise opti structure
opti = casadi.Opti();

% parameters
allparms = {'f','k11','k12','k21','k22','JF','koop','J1','J2', 'kon', 'koff', 'kse','kse0', 'kpe', 'Fpe0','b','k','dLcrit', 'gamma', 'kF', 'vmax', 'kappa', 'ps2', 'act_max'};

% create variables for all parameters
for i = 1:length(allparms)
    eval([allparms{i}, ' = ', num2str(parms.(allparms{i})),';'])
end

if ~isempty(optparms)
    % get bounds and normalized values
    lb = nan(1, length(optparms));
    ub = nan(1, length(optparms));
    nv = nan(1, length(optparms));

    for i = 1:length(optparms)
        lb(i) = bnds.(optparms{i})(1);
        ub(i) = bnds.(optparms{i})(2);
        nv(i) = (parms.(optparms{i}) - lb(i)) / (ub(i)-lb(i));
    end

    if (sum(nv > 1) + sum(nv < 0)) > 0
        disp('Warning: initial guess out of bounds')
        keyboard
    end

    % create opti variable for normalized parameters that are fitted
    s = opti.variable(1, length(optparms));

    % normalized values should be between 0 and 1
    opti.subject_to(0 < s < 1);

    % initial guess = current values
    opti.set_initial(s, nv);

    % recalculate non-normalized parameters from normalized parameters
    for i = 1:length(optparms)
        eval([optparms{i},' = s(', num2str(i), ')* ', num2str(ub(i)-lb(i)), '+', num2str(lb(i)), ';'])
    end

    JF = kF / J1;
end

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

%% define opti variables, specify constraints and initial guesses

if parms.f > 0 % biophysical models
    % define opti states (defined as above)
    Q0  = opti.variable(1,N);
    Q2  = opti.variable(1,N);
    L   = opti.variable(1,N);
%     Lts  = opti.variable(1,N);
    
    % define extra variables
    p   = opti.variable(1,N); % mean strain of the distribution
    q   = opti.variable(1,N); % standard deviation strain of the distribution
    Fse = opti.variable(1,N);
    Q1  = opti.variable(1,N);

    % (slack) controls (defined as above)
    dQ0dt  = opti.variable(1,N);
%     dQ1dt  = opti.variable(1,N);
    dQ2dt  = opti.variable(1,N);
    Ld     = opti.variable(1,N);
        
    % forces
    dlse = Lts - L;
    opti.subject_to(Fse == kse0 * (exp(kse*dlse)-1));
    
%     Fpe = kpe * log(1+exp(L*parms.K))/parms.K + Fpe0;
    Fpe = kpe * L + Fpe0;
    Fce = Fse - Fpe;
    
    opti.subject_to(Q1 == Fce - Q0);
    
%     opti.subject_to(Fce + Fpe == Fse);
 
    % instantaneous stiffnesses
%     Kpe = kpe .* (1 - 1./(exp(parms.K*L)+1));
    Kpe = kpe;
%     Kse = kse * (Fse + kse0);
    Kse = kse0 * kse * exp(kse*dlse);
    
    % extra constraints
    opti.subject_to(Q1 - Q0 .* p == 0);
    opti.subject_to(Q2 - Q0 .* (p.^2 + q) == 0);
    
    % (potentially) simple bounds
    opti.subject_to(q >= 0);
    opti.subject_to(Q0 > 0);
%     opti.subject_to(L < Lts); % avoid negative dlse
%     opti.subject_to(Q1 > -Q0);
%     opti.subject_to(Q2 > 0);
    opti.subject_to(Fse > 0);
%     opti.subject_to(Fce > 0);
    
    % initial guess states
    opti.set_initial(Q0, IG.Q0i);
    opti.set_initial(Q1, IG.Q1i);
    opti.set_initial(Q2, IG.Q2i);
    opti.set_initial(L, IG.Li);
%     opti.set_initial(Lts, IG.Lti);
      
    % initial guess extra variables
    opti.set_initial(p, IG.pi);
    opti.set_initial(q, IG.qi);
    opti.set_initial(Fse, IG.Fsei);
%     opti.set_initial(Fce, IG.Fcei);
    
    % initial guess slack controls
    opti.set_initial(dQ0dt, IG.dQ0dti);
%     opti.set_initial(dQ1dt, IG.dQ1dti);
    opti.set_initial(dQ2dt, IG.dQ2dti);
    opti.set_initial(Ld, IG.Ldi);

    % if thin filament cooperativity
    if parms.koop > 0
        Non = opti.variable(1,N);
        opti.subject_to(Non >= 0);
        opti.set_initial(Non, IG.Noni);
    end
    
    % if thick filament cooperativity
    if parms.J1 > 0    % cooperative models
        DRX = opti.variable(1,N);
        opti.subject_to(DRX >= 0);
        opti.set_initial(DRX, IG.DRXi);
    end
    
    % if forcibly detached
    if parms.b > 0 % FD model
        R = opti.variable(1,N);
        opti.subject_to(R >= 0);
        opti.set_initial(R, IG.Ri);
    else
        R = zeros(1, N);
    end
    
else % Hill-type model
    
    % CE force, length and velocity
    Fce  = opti.variable(1,N);
    L  = opti.variable(1,N);
    v  = opti.variable(1,N);
    
    opti.set_initial(Fce, IG.Fi);
    opti.set_initial(v, IG.vi);
    opti.set_initial(L, IG.Li);
    
    % activation from Ca
    Act = parms.act_max * Cas.^n ./ (kappa^n + Cas.^n);
    Act = log(1+exp(Act*parms.K))/parms.K; % avoid small numbers
    
    % activation-normalized force
    Frel = Fce ./ Act;
    
    % force-velocity relation
    vi = vmax/(2*parms.e(2))*(-exp(parms.e(4)/parms.e(1)-Frel/parms.e(1))+exp(Frel/parms.e(1)-parms.e(4)/parms.e(1))-2*parms.e(3));
    
    % interverted tendon stress-strain
    dLt = log(Fce./kse0 + 1) / kse;
    Li = Lts - dLt; % CE = FIBER - SE
    
end

%% dynamics constraints
if parms.f > 0 % biophysical models    

    
    if parms.koop == 0 % cooperative
        Non = Cas.^n ./ (kappa^n + Cas.^n); % sigmoid
    end
    
    if parms.J1 == 0
        DRX = 1 - Q0;
    end
    
    %% cross-bridge dynamics
    dQ1dt = 0;
    [error_Q0, ~, error_Q2, F0dot, dRdt] = MuscleEquilibrium(Q0, Q1, p, q, dQ0dt, dQ1dt, dQ2dt, f, parms.w, k11, k12, k21, k22, Non, Ld, DRX, b, k, R, dLcrit, ps2, parms.approx); 
    opti.subject_to(error_Q0(:) == 0);
%     opti.subject_to(error_Q1(:) == 0);
    opti.subject_to(error_Q2(:) == 0);
    
    % needed for slack controls
    opti.subject_to((dQ0dt(1:N-1) + dQ0dt(2:N))*dt/2 + Q0(1:N-1) == Q0(2:N));
%     opti.subject_to((dQ1dt(1:N-1) + dQ1dt(2:N))*dt/2 + Q1(1:N-1) == Q1(2:N));
    opti.subject_to((dQ2dt(1:N-1) + dQ2dt(2:N))*dt/2 + Q2(1:N-1) == Q2(2:N));
    
    %% length dynamics
    error_length    = LengthEquilibrium(Q0, F0dot, Ld, vts, Kse, Kpe, parms.gamma);
    opti.subject_to(error_length(:) == 0);
    
%     V = vts * parms.gamma;
    
    % needed for slack controls
    opti.subject_to((Ld(1:N-1) + Ld(2:N))*dt/2 + L(1:N-1) == L(2:N));
%     opti.subject_to((V(1:N-1) + V(2:N))*dt/2 + Lts(1:N-1) == Lts(2:N));
%     opti.subject_to(Lts(1) == 0);
    
    %% cooperativity
    if parms.koop > 0 % cooperative
        [Jon, Joff] = ThinFilament_Dynamics(Cas, Q0, Non, kon, koff, koop, parms.Noverlap);
        dNondt = Jon - Joff;
        opti.subject_to((dNondt(1:N-1) + dNondt(2:N))*dt/2 + Non(1:N-1) == Non(2:N));
    end
    
    if parms.J1 > 0
        [k1, k2] = ThickFilament_Dynamics(Q0, Fce, DRX, J1, J2, JF, parms.Noverlap, R);
        dDRXdt = k1 - k2 - (dQ0dt + dRdt);
        opti.subject_to((dDRXdt(1:N-1) + dDRXdt(2:N))*dt/2 + DRX(1:N-1) == DRX(2:N));
    end
    
   %% forcibly detached states
    if parms.b > 0
%         opti.subject_to(error_R(:) == 0);
        opti.subject_to((dRdt(1:N-1) + dRdt(2:N))*dt/2 + R(1:N-1) == R(2:N));
    end

else % Hill-type
    
    % error terms
    opti.subject_to(vi - v == 0);
    opti.subject_to(Li - L == 0);
    opti.subject_to((v(1:N-1) + v(2:N))*dt/2 + L(1:N-1) == L(2:N));
end

%% cost
Frel = Fse * parms.Fscale;

% not needed for Hill-type, because already enforced by dynamics
if parms.f > 0 % biophysical models
    opti.subject_to(Frel(1:10) == 1);
end

Fcost = (Frel(idF) - Fts(idF)).^2;

% cost function
J = 0;
J = J + w(1) * sum(Fcost); % force-velocity fitting

if parms.f > 0
    J = J + w(3) * (sum(dQ0dt(idC).^2) + sum(dQ2dt(idC).^2)); % regularization term
else
    J = J + w(3) * (sum(vi(idC).^2)); % regularization term
end

% optimize
opti.minimize(J);

%% Solve problem
% options for IPOPT

options.ipopt.linear_solver = 'mumps';
% options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.mu_strategy           = 'adaptive';
options.detect_simple_bounds           = true;

options.ipopt.max_iter           = 500;

opti.solver('ipopt',options);
%%

% error_length    = LengthEquilibrium(Q0, F0dot, Ld, vts, Kse, Kpe, parms.gamma);
% [error_Q0, ~, error_Q2, F0dot, dRdt] = MuscleEquilibrium(Q0, Q1, p, q, dQ0dt, dQ1dt, dQ2dt, f, parms.w, k11, k12, k21, k22, Non, Ld, DRX, b, k, R, dLcrit, ps2, parms.approx); 
%     
% x = opti.debug.value(Frel,opti.initial())
% y = opti.debug.value(error_Q2,opti.initial())
% 
% close all
% figure(1)
% subplot(211)
% plot(x)
% 
% subplot(212)
% plot(y); hold on
% % plot(IG.F0doti,'--')
%%

% Solve the OCP
% p_opts = struct('detect_simple_bounds', true);
% s_opts = struct('max_iter', 500);
% opti.solver('ipopt',p_opts,s_opts);

% visualize
% opti.callback(@(i) plot(toc, [Fts; opti.debug.value(Frel)]))

try
    sol = opti.solve();
    
    out.F     = sol.value(Frel);
    out.J     = sol.value(J);
    out.Fcost = sol.value(Fcost);
    out.s = sol.value(s);
    
    if parms.f > 0
        % Extract the result
        out.Q0    = sol.value(Q0);
        out.Q1    = sol.value(Q1);
        out.Q2    = sol.value(Q2);
        out.Non    = sol.value(Non);
        out.DRX    = sol.value(DRX);
        
        out.dQ0dt = sol.value(dQ0dt);
        out.dQ1dt = sol.value(dQ1dt);
        out.dQ2dt = sol.value(dQ2dt);
        
        out.Lts = sol.value(Lts);
        
        out.dDRXdt = sol.value(dDRXdt);
        out.dNondt = sol.value(dNondt);
        
        out.F0dot = sol.value(F0dot);
        out.p = sol.value(p);
        out.q = sol.value(q);
        out.L = sol.value(L);
        out.Ld = sol.value(Ld);
        out.Fce = sol.value(Fce);
        out.Fse = sol.value(Fse);
        out.Fpe = sol.value(Fpe);
    end
    
    % extract the parameters
    for i = 1:length(optparms)
        parms.(optparms{i}) = eval(['sol.value(',optparms{i},');']);
    end
    
catch
    sol = opti.debug();
    
    out.F     = opti.debug.value(Frel);
    out.J     = opti.debug.value(J);
    out.Fcost = opti.debug.value(Fcost);
    
    if parms.f > 0
        % Extract the result
        out.Q0    = opti.debug.value(Q0);
        out.Q1    = opti.debug.value(Q1);
        out.Q2    = opti.debug.value(Q2);
        out.dQ0dt = opti.debug.value(dQ0dt);
        out.dQ1dt = opti.debug.value(dQ1dt);
        out.dQ2dt = opti.debug.value(dQ2dt);
    end
    
    % extract the parameters
    for i = 1:length(optparms)
        parms.(optparms{i}) = eval(['opti.debug.value(',optparms{i},');']);
    end
end

if parms.f > 0
    out.Fdot  = out.dQ0dt + out.dQ1dt;
end

out.t     = 0:dt:(N-1)*dt;

end
