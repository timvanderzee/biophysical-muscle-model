function[parms, out, opti] = fit_model_parameters_v2(optparms, w, data, parms, IG, bnds)

import casadi.*

% initialise opti structure
opti = casadi.Opti();

% parameters
allparms = {'f','k11','k12','k21','k22','JF','koop','J1','J2', 'kon', 'koff', 'kse','kse0', 'kpe', 'Fpe0','b','k','dLcrit', 'gamma', 'kF', 'vmax', 'kappa', 'ps2', 'act_max', 'Lce0'};

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
% define opti states (defined as above)
Q0  = opti.variable(1,N); % derivative constraint
Q2  = opti.variable(1,N); % derivative constraint
% Fse = opti.variable(1,N); % derivative constraint
Q1  = opti.variable(1,N);
Non = opti.variable(1,N);  % derivative constraint
DRX = opti.variable(1,N);  % derivative constraint

Ld = opti.variable(1,N);  % derivative constraint
opti.set_initial(Ld, IG.Ldi);

% p and q
p   = opti.variable(1,N); % mean strain of the distribution
q   = opti.variable(1,N); % standard deviation strain of the distribution
% 
% y   = opti.variable(1); % standard deviation strain of the distribution
% opti.set_initial(y, 1);

% compute other things from force
% dlse = log(Fse/kse0+1)/kse;
% L = Lts - dlse;
% Fpe = kpe * L + Fpe0;
Fpe = 0;

Fce = Q0 + Q1;
Fse = Fce;
Kse = kse * (Fse + kse0);
% Kpe = kpe;
% Kpe = 0;

% extra constraints
% opti.subject_to(Q1 + Q0 + Fpe == Fse);
opti.subject_to(Q1 - Q0 .* p == 0);
opti.subject_to(Q2 - Q0 .* (p.^2 + q) == 0);

% bounds
opti.subject_to(q > 0);
opti.subject_to(p > -5);
opti.subject_to(p < 5);
opti.subject_to(Q0 > 0);
opti.subject_to(Fse > 0);

% initial guess states
opti.set_initial(Q0, IG.Q0i);
opti.set_initial(Q1, IG.Q1i);
opti.set_initial(Q2, IG.Q2i);

% initial guess extra variables
opti.set_initial(p, IG.pi);
opti.set_initial(q, IG.qi);
% opti.set_initial(Fse, IG.Fsei);

opti.subject_to(Non > 0);
opti.set_initial(Non, IG.Noni);

opti.subject_to(DRX > 0);
opti.set_initial(DRX, IG.DRXi);

R = 0;
dRdt = 0;

%% Passive force
% low risk: only add Fpe at the end
dlse = log(Fse/kse0+1)/kse;
L = Lts - dlse;
dLce = L - Lce0;
K = 1;
Lcor = log(1+exp(dLce*K))/K;



% passive force
Fp = kpe * Lcor;

%% cross-bridge dynamics
% points where integrals is evaluated
k1 = [k11 k12];
k2 = [k21 -k22];

% get functions
[IG, IGef] = get_IG_IGEf(parms.approx);

% Compute Qdot
[Q0dot, Q1dot, Q2dot] = CrossBridge_Dynamics(Q0, p, q, f, parms.w, k1, k2, IGef, Non, DRX, IG, b, k, R, dLcrit, ps2);

% velocity - independent derivative
F0dot  = Q1dot + Q0dot;

% velocity
opti.subject_to(Ld .* (Q0 + Kse + Kpe) == (vts .* parms.gamma .* Kse - F0dot));

% cross-bridge derivatives
dQ0dt = Q0dot;
dQ1dt = Q1dot + 1 * Ld .* Q0;
dQ2dt = Q2dot + 2 * Ld .* Q1;

% force-rate
% dFcedt = dQ1dt + dQ0dt;
% dFpedt = Ld * Kpe;
% dFsedt = dFcedt + dFpedt;

% enforce dynamic constraints
opti.subject_to((dQ0dt(1:N-1) + dQ0dt(2:N))*dt/2 + Q0(1:N-1) == Q0(2:N));
opti.subject_to((dQ1dt(1:N-1) + dQ1dt(2:N))*dt/2 + Q1(1:N-1) == Q1(2:N));
opti.subject_to((dQ2dt(1:N-1) + dQ2dt(2:N))*dt/2 + Q2(1:N-1) == Q2(2:N));
% opti.subject_to((dFsedt(1:N-1) + dFsedt(2:N))*dt/2 + Fse(1:N-1) == Fse(2:N));

%% cooperativity
[Jon, Joff] = ThinFilament_Dynamics(Cas, Q0, Non, kon, koff, koop, parms.Noverlap);
dNondt = Jon - Joff;
opti.subject_to((dNondt(1:N-1) + dNondt(2:N))*dt/2 + Non(1:N-1) == Non(2:N));

[k1, k2] = ThickFilament_Dynamics(Q0, Fce, DRX, J1, J2, JF, parms.Noverlap, R);
dDRXdt = k1 - k2 - (dQ0dt + dRdt);
opti.subject_to((dDRXdt(1:N-1) + dDRXdt(2:N))*dt/2 + DRX(1:N-1) == DRX(2:N));

%% cost
Frel = (Fse + Fp) * parms.Fscale;

% not needed for Hill-type, because already enforced by dynamics
opti.subject_to(Frel(1) == 1);

Fcost = (Frel(idF) - Fts(idF)).^2;

% cost function
J = 0;
J = J + w(1) * sum(Fcost); % force-velocity fitting
J = J + w(3) * (sum(dQ0dt(idC).^2) + sum(dQ1dt(idC).^2) + sum(dQ2dt(idC).^2)); % regularization term

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

try
    sol = opti.solve();
    
    out.F     = sol.value(Frel);
    out.J     = sol.value(J);
    out.Fcost = sol.value(Fcost);
    
    if exist('s')
        out.s = sol.value(s);
    end
    
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
        %         out.L = sol.value(L);
        %         out.Ld = sol.value(Ld);
        out.Fce = sol.value(Fce);
        out.Fse = sol.value(Fse);
        out.Fpe = sol.value(Fpe);
    end
    
    % extract the parameters
    for i = 1:length(optparms)
        parms.(optparms{i}) = eval(['sol.value(',optparms{i},');']);
    end
    
catch
    disp('Optimal solution not found')
    sol = opti.debug();
    
    out.F     = opti.debug.value(Frel);
    out.J     = opti.debug.value(J);
    out.Fcost = opti.debug.value(Fcost);
    out.Lts = sol.value(Lts);
    
    if parms.f > 0
        % Extract the result
        out.Q0    = opti.debug.value(Q0);
        out.Q1    = opti.debug.value(Q1);
        out.Q2    = opti.debug.value(Q2);
        
        out.Non    = opti.debug.value(Non);
        out.DRX    = opti.debug.value(DRX);
        
        out.dQ0dt = opti.debug.value(dQ0dt);
        out.dQ1dt = opti.debug.value(dQ1dt);
        out.dQ2dt = opti.debug.value(dQ2dt);
    end
    
    % extract the parameters
    for i = 1:length(optparms)
        parms.(optparms{i}) = eval(['opti.debug.value(',optparms{i},');']);
    end
end

out.t     = 0:dt:(N-1)*dt;

end
