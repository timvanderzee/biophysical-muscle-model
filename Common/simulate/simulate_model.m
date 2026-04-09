function[sol, out] = simulate_model(model, odefunc, modelfunc, input, newparms)
% out = struct();

if strcmp(model, 'Hill-type no SE')
    % CE activation
    for i = 1:length(input) % loop over phases
        newparms.actfunc = @(Ca,parms)parms.act_max*Ca.^parms.n./(parms.kappa^parms.n+Ca.^parms.n);
        a = newparms.actfunc(input(i).Cas, newparms);
        vce = input(i).v * newparms.gamma;
        
        % PE
        dLce = input(i).L * newparms.gamma - newparms.Lce0;
        Fpe = newparms.kpe * dLce;
        
        % force-velocity
        Fce = a .* newparms.vF_func(vce, newparms);
        
        out(i).F = (Fpe + Fce) * newparms.Fscale;
        out(i).t = input(i).t;
        sol = [];
    end
    
else
    
    % determine initial state
    % we need to define the model states at t = 0
    if ~isfield(newparms, 'X0')
        [X0, newparms] = get_steady_state(model, odefunc, modelfunc, newparms, input(1).Cas(1),input(1).L(1)*newparms.gamma);
    else
        X0 = newparms.X0;
        newparms.vts = [0 0];
        newparms.ti = [0 10];
        newparms.Lts = input(1).Cas(1) * [1 1];
        newparms.Cas = input(1).L(1)*newparms.gamma * [1 1];

    end
    
    if contains(char(odefunc), 'i') % implicit
        sol = odefunc(@(t,y,yp) modelfunc(t,y,yp, newparms), [0 .001], X0, zeros(size(X0)), []);
    else
        sol = odefunc(@(t,y) modelfunc(t,y, newparms), [0 .001], X0, []); % just to define sol
    end
    
    % simulate model
    for i = 1:length(input) % loop over phases
        
        % obtain the inputs for this phase
        newparms.ti     = input(i).t;
        newparms.Lts    = input(i).L * newparms.gamma;
        newparms.vts    = input(i).v;
        newparms.Cas    = input(i).Cas;
        
        % simulate
        if contains(char(odefunc), 'i') % implicit
            sol(i) = odefunc(@(t,y,yp) modelfunc(t,y,yp, newparms), [0 max(newparms.ti)], X0, zeros(size(X0)), odeset('maxstep', 1e-2));
        else
            sol(i) = odefunc(@(t,y,yp) modelfunc(t,y, newparms), [0 max(newparms.ti)], X0, odeset('maxstep', 1e-2));
        end
        
        % save final state (used as initial state for next simulation)
        X0 = sol(i).y(:,end);
    end
    
    % get force
    for i = 1:length(input) % loop over phases
        out(i) = get_forces_from_state(model, modelfunc, sol(i), input(i), newparms);
    end
end
end

function[out] = get_forces_from_state(model, modelfunc, sol, input, newparms)

t   = sol.x;

if contains(model, 'XB')
    
    if contains(char(modelfunc), 'discretized')
        Lce = sol.y(end-3,:);
        
        dlse = interp1(input.t, input.L * newparms.gamma, t) - Lce;
        F = newparms.Fse_func(dlse, newparms);

        FXB = nan(1, length(sol.x));
        
        for iiii = 1:length(sol.x)
            n = sol.y(1:end-4,iiii);
            xi = newparms.xi + Lce(iiii);
            FXB(iiii) = trapz(xi, xi .* n') + trapz(xi, n');
        end
        
        [Fpe, ~] = get_PE(Lce, newparms);

        Fse = (FXB + Fpe) * newparms.Fscale;
                                
    else
        
        Fse = sol.y(3,:) * newparms.Fscale;
        dLse = max(newparms.Lse_func(Fse/newparms.Fscale, newparms), 0); % can't be negative
        Lce = interp1(input.t, input.L * newparms.gamma, t) - dLse;
    end
else
    Lce = sol.y;
    dLse = interp1(input.t, input.L * newparms.gamma, t) - Lce;
    Fse = newparms.Fse_func(dLse, newparms) * newparms.Fscale;
end

if ~isfield(newparms, 'Lce0')
    newparms.Lce0 = 0;
end

dLce = Lce - newparms.Lce0;
Fpe = newparms.kpe * dLce  * newparms.Fscale;
Fpe(newparms.K*dLce < 10) = newparms.kpe * log(1+exp(dLce(newparms.K*dLce < 10)*newparms.K))/newparms.K  * newparms.Fscale;

Fce = Fse - Fpe;

out.t = t;
out.F = Fse;
out.Fpe = Fpe;
out.Fce = Fce;
out.Lce = Lce;
out.dLse = dLse;
% out(i).t, out(i).F, out(i).Fpe, out(i).Fce, out(i).Lce
end


