function[parms, out, opti] = fit_model_parameters_passive(optparms, w, data, parms, IG, bnds)

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
   
end

%% extract input and target
Fts = data.F;
toc = data.t;
Lts = data.L;
idF = data.idF;

N = length(toc);
dt = mean(diff(toc));

Fse = opti.variable(1,N); % derivative constraint
opti.set_initial(Fse, IG.F);
opti.subject_to(Fse > 0);

% lengths
dlse = log(Fse/kse0+1)/kse;
L = Lts - dlse;
dLce = L - Lce0;

% parallel force
Fpe = kpe * dLce;

% constraint (assume Fce = 0)
opti.subject_to(Fse == Fpe);

% cost
Frel = Fse * parms.Fscale;
Fcost = (Frel(idF) - Fts(idF)).^2;

% cost function
J = 0;
J = J + sum(Fcost); % force-velocity fitting

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

%% get outputs

try
    sol = opti.solve();
    
    out.F     = sol.value(Frel);
    out.J     = sol.value(J);
    out.Fcost = sol.value(Fcost);
    
    if exist('s')
        out.s = sol.value(s);
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
    
        if exist('s')
        out.s = sol.value(s);
        end
    
    % extract the parameters
    for i = 1:length(optparms)
        parms.(optparms{i}) = eval(['opti.debug.value(',optparms{i},');']);
    end
end

out.t     = 0:dt:(N-1)*dt;

end
