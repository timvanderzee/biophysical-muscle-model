function[sol, out] = simulate_model(model, modelfunc, input, newparms)
    out = struct();
    
   if strcmp(model, 'Hill-type no SE')
    % CE activation
      for i = 1:length(input) % loop over phases
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
    X0 = get_steady_state(model, modelfunc, newparms, input(1).Cas(1));
    sol = ode15s(@(t,y,yp) modelfunc(t,y, newparms), [0 .001], X0, []); % just to define sol

    % simulate model
    for i = 1:length(input) % loop over phases

        % obtain the inputs for this phase
        newparms.ti     = input(i).t;
        newparms.Lts    = input(i).L * newparms.gamma;
        newparms.vts    = input(i).v;
        newparms.Cas    = input(i).Cas;

        % simulate
        sol(i) = ode15s(@(t,y,yp) modelfunc(t,y, newparms), [0 max(newparms.ti)], X0, odeset('maxstep', 1e-2));

        % save final state (used as initial state for next simulation)
        X0 = sol(i).y(:,end);
    end

    % get force
    for i = 1:length(input) % loop over phases
        [out(i).t, out(i).F, ~, ~] = get_forces_from_state(model, sol(i), input(i), newparms);
    end
   end
end

function[t, Fse, Fpe, Fce] = get_forces_from_state(model, sol, input, newparms)

    t   = sol.x;
    
    if contains(model, 'XB')
        Fse = sol.y(3,:) * newparms.Fscale;
        dLse = max(newparms.Lse_func(Fse, newparms), 0); % can't be negative
        Lce = interp1(input.t, input.L * newparms.gamma, t) - dLse;
    else
        Lce = sol.y;
        dLse = interp1(input.t, input.L * newparms.gamma, t) - Lce;
        Fse = newparms.Fse_func(dLse, newparms) * newparms.Fscale;
    end
    
    dLce = Lce - newparms.Lce0;
    Fpe = newparms.kpe * dLce  * newparms.Fscale;
    Fpe(newparms.K*dLce < 10) = newparms.kpe * log(1+exp(dLce(newparms.K*dLce < 10)*newparms.K))/newparms.K  * newparms.Fscale;

    Fce = Fse - Fpe;
end


