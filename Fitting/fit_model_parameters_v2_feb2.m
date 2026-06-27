function[parms, out, opti] = fit_model_parameters_v2_feb2(model, parms, optparms, bnds, data, IG, weights)

%% extract input and target
% mandatory
Cas = data.Cas;                 % calcium concentration
vts = data.v;                   % crossbridge velocity
Fts = data.F;                   % fiber force
toc = data.t;                   % time vector
Lts = data.L * parms.gamma;     % crossbridge length change

% optional: active versus passive indices
if isfield(data, 'idA'), idA = data.idA;
else,                    idA = 1:length(Cas); % assume everything is active
end

if isfield(data, 'idP'), idP = data.idP;
else,                    idP = []; % assume nothing is passive
end

% optional: cost function indices
if isfield(data, 'idF'), idF  = data.idF; % active force term
else,                    idF = find(isfinite(Fts(idA)));
end

if isfield(data, 'idFP'), idFP = data.idFP; % passive force term
else,                     idFP = idP;
end

if isfield(data, 'idC'),  idC  = data.idC; % regularization term
else,                     idC = [];
end

N = length(idA);
dt = mean(diff(toc));
out.t     = 0:dt:(N-1)*dt;

%% normalize parameters
import casadi.*
opti = casadi.Opti(); % initialise opti structure

% create variables for all parameters
allparms = fieldnames(parms);
for i = 1:length(allparms)
    if isscalar(parms.(allparms{i})) && isnumeric(parms.(allparms{i}))
        eval([allparms{i}, ' = ', num2str(parms.(allparms{i})),';'])
    end
end

% create variable for parameters that will be fitted
if ~isempty(optparms)
    % get bounds and normalized values
    lb = nan(1, length(optparms));
    ub = nan(1, length(optparms));
    nv = nan(1, length(optparms));
    
    % convert into normalized value between bounds
    for i = 1:length(optparms)
        lb(i) = bnds.(optparms{i})(1);
        ub(i) = bnds.(optparms{i})(2);
        nv(i) = (parms.(optparms{i}) - lb(i)) / (ub(i)-lb(i));
    end
    
    % make sure that initial guess is not at the bounds
    if (sum(nv > 1) + sum(nv < 0)) > 0
        disp('Warning: initial guess out of bounds. Adjusting ...')
        nv(nv > 1) = .9;
        nv(nv < 0) = .1;
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

%% passive force fitting (pCa 9)
if ~isempty(idP)
    NP = length(idP);
    
    FseP = opti.variable(1,NP);
    opti.set_initial(FseP, IG.F9(idP));
    opti.subject_to(FseP > 0);
    
    % lengths
    dlseP = log(FseP/kse0+1)/kse;
    LceP  = Lts(idP) - dlseP;
    dLceP = LceP - Lce0;
    
    % parallel force
    FpeP = kpe * dLceP;
    
    % constraint (assume Fce = 0)
    opti.subject_to(FseP == FpeP);
    
    % cost
    FrelP = FseP * parms.Fscale;
    FcostP = (FrelP(idFP - N) - Fts(idFP)).^2;
end

%% active force fitting: define opti variables, specify constraints and initial guesses
if contains(model, 'XB')
    % define opti variables
    Q0  = opti.variable(1,N); % XBs attached
    Q1  = opti.variable(1,N); % XB force
    Q2  = opti.variable(1,N); % XB elastic strain energy
    Fse = opti.variable(1,N); % SE force
    p   = opti.variable(1,N); % mean strain of the distribution
    q   = opti.variable(1,N); % standard deviation strain of the distribution
    Ld  = opti.variable(1,N);  % % fiber velocity
    
    % constraints that specify relation between opti variables
    opti.subject_to(Q1 - Q0 .* p == 0);
    opti.subject_to(Q2 - Q0 .* (p.^2 + q) == 0);
    
    % bounds
    opti.subject_to(q > 0);
    opti.subject_to(p > -5);
    opti.subject_to(p < 5);
    opti.subject_to(Q0 > 1e-6);
    opti.subject_to(Fse > 0);
    % opti.subject_to(Q0+Q1 > 0);
    
    % initial guess
    opti.set_initial(Q0, IG.Q0i(idA));
    opti.set_initial(Q1, IG.Q1i(idA));
    opti.set_initial(Q2, IG.Q2i(idA));
    opti.set_initial(Ld, IG.Ldi(idA));
    opti.set_initial(p, IG.pi(idA));
    opti.set_initial(q, IG.qi(idA));
    opti.set_initial(Fse, IG.Fsei(idA));
    
    %% active force fitting: forces and stiffnesses
    % passive force
    dLse = log(Fse/kse0+1)/kse;
    Lce = Lts(idA) - dLse;
    dLce = Lce - Lce0;
    K = 1;
    Fpe = kpe * log(1+exp(dLce*K))/K;
    
    % active force
    Fce = Q0 + Q1;
    
    % stiffneses
    Kse = kse * (Fse + kse0);
    Kpe = kpe .* (1 - 1./(exp(K*dLce)+1));
    
    % force constraint
    opti.subject_to(Fse == Fce + Fpe);
    
    % scale the force
    Frel = Fse * parms.Fscale;
    
    % optional: constrain initial value of the force
    % opti.subject_to(Frel(1) == 1);
    opti.subject_to(Frel(1) == Fts(idF(1)));
    
    % force-term in cost function
    Fcost = (Frel(idF) - Fts(idF)).^2;
    
    %% active force fitting: optional additional variables
    % default: assume regular 2-state model
    Non = Cas(idA).^n ./ (kappa^n + Cas(idA).^n);
    DRX = 1 - Q0;
    R   = 0;
    
    % thin filament activation
    if contains(model, 'coop')
        Non = opti.variable(1,N);  % derivative constraint
        opti.subject_to(Non > 0);
        opti.subject_to(Non < 2);
        opti.set_initial(Non, IG.Noni(idA));
    end
    
    % thick filament activation
    if contains(model, '3-state') || contains(model, '4-state')
        DRX = opti.variable(1,N);  % derivative constraint
        opti.subject_to(DRX > 0);
        opti.subject_to(DRX < 2);
        opti.set_initial(DRX, IG.DRXi(idA));
    end
    
    % forcibly detached state
    if contains(model, '4-state')
        R = opti.variable(1,N);  % derivative constraint
        %    opti.subject_to(R > 0);
        opti.set_initial(R, IG.Ri(idA));
    end
    
    %% active force fitting: cross-bridge dynamics
    % points where integrals is evaluated
    k1 = [k11 k12];
    k2 = [k21 -k22];
    
    % get functions
    [IG, IGef] = get_IG_IGEf(parms.approx);
    
    % Compute Qdot
    [Q0dot, Q1dot, Q2dot, Rdot] = CrossBridge_Dynamics(Q0, p, q, f, parms.w, k1, k2, IGef, Non, DRX, IG, b, k, R, dLcrit, ps2);
    
    % velocity - independent derivative
    F0dot  = Q1dot + Q0dot;
    
    % velocity constraint
    opti.subject_to(Ld .* (Q0 + Kse + Kpe) == (vts(idA) .* parms.gamma .* Kse - F0dot));
    
    % cross-bridge derivatives
    dQ0dt = Q0dot;
    dQ1dt = Q1dot + 1 * Ld .* Q0;
    dQ2dt = Q2dot + 2 * Ld .* Q1;
    
    % enforce dynamic constraints
    opti.subject_to((dQ0dt(1:N-1) + dQ0dt(2:N))*dt/2 + Q0(1:N-1) == Q0(2:N));
    opti.subject_to((dQ1dt(1:N-1) + dQ1dt(2:N))*dt/2 + Q1(1:N-1) == Q1(2:N));
    opti.subject_to((dQ2dt(1:N-1) + dQ2dt(2:N))*dt/2 + Q2(1:N-1) == Q2(2:N));
    
    %% active force fitting: additional dynamics
    % forcibly detached state
    if contains(model, '4-state')
        opti.subject_to((Rdot(1:N-1) + Rdot(2:N))*dt/2 + R(1:N-1) == R(2:N));
    end
    
    % thin filament dynamics
    if contains(model, 'coop')
        [Jon, Joff] = ThinFilament_Dynamics(Cas(idA), Q0, Non, kon, koff, koop, parms.Noverlap);
        dNondt = Jon - Joff;
        opti.subject_to((dNondt(1:N-1) + dNondt(2:N))*dt/2 + Non(1:N-1) == Non(2:N));
    end
    
    % thick filament dynamics
    if contains(model, '3-state') || contains(model, '4-state')
        [k1, k2] = ThickFilament_Dynamics(Q0, Fce, DRX, J1, J2, JF, parms.Noverlap, R);
        dDRXdt = k1 - k2 - (dQ0dt + Rdot);
        opti.subject_to((dDRXdt(1:N-1) + dDRXdt(2:N))*dt/2 + DRX(1:N-1) == DRX(2:N));
    end
    
    
else % Hill-type
    
    % variables
    Fse = opti.variable(1, N);
    opti.set_initial(Fse, IG.Fsei(idA));
    
    % CE activations
    a = act_max * Cas(idA).^n ./ (kappa^n + Cas(idA).^n);
    
    if strcmp(model, 'Hill-type SE')
        vce = opti.variable(1, N);
        opti.set_initial(vce, IG.vcei(idA));

        % PE
        dLse = log(Fse/kse0+1)/kse;
        Lce = Lts(idA) - dLse;
        dLce = Lce - Lce0;
        K = 1;
        Fpe = kpe * log(1+exp(dLce*K))/K;

        % force-velocity
        Fce = a .* (parms.e(1)*log((parms.e(2)*vce./vmax+parms.e(3))+sqrt((parms.e(2)*vce./vmax+parms.e(3)).^2+1))+parms.e(4));
       
        % velocity is time derivative of length
        opti.subject_to((vce(1:N-1) + vce(2:N))*dt/2 + Lce(1:N-1) == Lce(2:N)); 
    
    else % no SE
        
        vce = vts(idA) * parms.gamma;
        
        % PE
        dLce = Lts(idA) - Lce0;
        K = 1;
        Fpe = kpe * log(1+exp(dLce*K))/K;

        % force-velocity
        Fce = a .* (parms.e(1)*log((parms.e(2)*vce./vmax+parms.e(3))+sqrt((parms.e(2)*vce./vmax+parms.e(3)).^2+1))+parms.e(4));
        
    end
    
    % force constraint
    opti.subject_to(Fce + Fpe == Fse);
   
    % Total force
    Frel = Fse * parms.Fscale;
    
    opti.subject_to(Frel(1) == Fts(idF(1)));
    
    % force-term in cost function
    Fcost = (Frel(idF) - Fts(idF)).^2;
end



%% specify cost to be minimized
% cost function
J = 0;
J = J + weights(1) * sum(Fcost); % active force cost

if ~isempty(idP) && length(weights)>1
    J = J + weights(2) * sum(FcostP);  % passive force cost
end

if ~isempty(idC) && length(weights)>2
    J = J + weights(3) * (sum(dQ0dt(idC).^2) + sum(dQ1dt(idC).^2) + sum(dQ2dt(idC).^2)); % regularization term
end

% minimize this cost function
opti.minimize(J);

%% optimization settings
options.detect_simple_bounds    = true;

% options for IPOPT
options.ipopt.linear_solver     = 'mumps';
options.ipopt.mu_strategy       = 'adaptive';
options.ipopt.max_iter          = 1e3;
opti.solver('ipopt',options);

%% solve the problem
try
    sol = opti.solve();
catch
    sol = opti.debug();
end

out.F       = sol.value(Frel);
if contains(model, 'XB')
    %% obtain output
    % extract variables
    out.Q0      = sol.value(Q0);
    out.Q1      = sol.value(Q1);
    out.Q2      = sol.value(Q2);
    out.Non     = sol.value(Non);
    out.DRX     = sol.value(DRX);
    out.dQ0dt   = sol.value(dQ0dt);
    out.dQ1dt   = sol.value(dQ1dt);
    out.dQ2dt   = sol.value(dQ2dt);
    out.Lts     = sol.value(Lts);
    out.F0dot   = sol.value(F0dot);
    out.p       = sol.value(p);
    out.q       = sol.value(q);
    out.Ld      = sol.value(Ld);
    out.Fce     = sol.value(Fce);
    out.Fse     = sol.value(Fse);
    out.Fpe     = sol.value(Fpe);
    out.s       = sol.value(s);
    out.J       = sol.value(J);
    out.Fcost   = sol.value(Fcost);
    
    if ~isempty(idP)
        out.Fcost9  = sol.value(FcostP);
    end
    
    % forcible detached
    if contains(model, '4-state')
        out.R = sol.value(R);
    end
    
    % thick filament cooperativity
    if contains(model, '3-state') || contains(model, '4-state')
        out.dDRXdt = sol.value(dDRXdt);
    end
    
    % thin filament cooperativity
    if contains(model, 'coop')
        out.dNondt = sol.value(dNondt);
    end

else % Hill-type specific
    out.a = sol.value(a);
end

% extract the parameters
for i = 1:length(optparms)
    parms.(optparms{i}) = eval(['sol.value(',optparms{i},');']);
end
  
% make sure that JF is corrected
if parms.J1 > 0
    parms.JF = parms.kF / parms.J1;
end

% time vector
out.t     = 0:dt:(N-1)*dt;

end
