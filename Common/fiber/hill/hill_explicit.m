function[xdot, Fce] = hill_explicit(t, x, parms)

    % states
    Lce = x(1);
   
    if numel(parms.Cas) == 1
        Ca = parms.Cas;
    else
        Ca = interp1(parms.ti, parms.Cas, t);
    end

    if numel(parms.Lts) == 1
        Lts = parms.Lts;
    else
        Lts = interp1(parms.ti, parms.Lts, t);
    end

    % activation from Ca
    parms.actfunc = @(Ca,parms)parms.act_max*Ca.^parms.n./(parms.kappa^parms.n+Ca.^parms.n);
    a = parms.actfunc(Ca, parms);
    a(a<parms.amin) = parms.amin;

    % elastic elements
    Lse = Lts - Lce;    

    Fpe = 0;
    
    Fse = parms.Fse_func(Lse, parms);
    Fse(Lse < 0) = 0;

    % contractile element
    Fce = Fse - Fpe;
    Fce(Fce < 0) = 0;

    % force-velocity
    Frel = Fce ./ a;
    Vce = parms.Fv_func(Frel, parms);

    % state derivative
    xdot = Vce;

end