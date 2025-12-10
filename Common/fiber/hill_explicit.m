function[xdot, Fce] = hill_explicit(t, x, parms, Ca)

    % states
    lce = x(1);
%     lmtc = x(2);

    lmtc = parms.Lts;

    % activation from Ca
    a = parms.actfunc(Ca, parms);
    a(a<parms.amin) = parms.amin;

    % elastic elements
    Lse = lmtc - lce;    
% 
%     if isfield(parms, 'Fpece_func')
%         Fpe = parms.Fpece_func(lce, parms);
%     else
%         Fpe = 0;
%     end

    Fpe = 0;
    
    Fse = parms.Fse_func(Lse, parms);
    Fse(Lse < 0) = 0;

    % contractile element
    Fce = Fse - Fpe;
    Fce(Fce < 0) = 0;

    % force-velocity
    Frel = Fce ./ a;
    vce = parms.Fv_func(Frel, parms);

    % state derivative
    xdot = [vce];

end