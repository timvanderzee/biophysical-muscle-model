function[F] = get_force(x, parms)

    L = x(end-3,:);

    Fpe = parms.Fpe_func(L, parms);
    
    
    parms.Fse_func = @(dlse, parms) parms.kse0*(exp(parms.kse*dlse)-1);
    F = parms.Fse_func(dlse, parms);
    
end