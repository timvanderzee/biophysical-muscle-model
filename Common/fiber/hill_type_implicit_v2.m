function[error, Fse] = hill_type_implicit_v2(t, x, xdot, parms)

% gamma = (.5*parms.s)/parms.h;

if numel(parms.Cas) == 1
    Ca = parms.Cas;
else
    Ca      = interp1(parms.ti, parms.Cas, t);
end

lmtc    = interp1(parms.ti, parms.Lts, t) * parms.gamma;

% states
lce = x(1);
vce = xdot(1);

% activation from Ca
a = parms.actfunc(Ca, parms);
a(a<parms.amin) = parms.amin;

% elastic elements
Lse = lmtc - lce;
Fse = parms.Fse_func(Lse, parms);
Fse(Lse < 0) = 0;

% activation-normalized force
Fce_rel = Fse ./ a;

vi = parms.Fv_func(Fce_rel, parms);

% error terms
error(1,1) = vce - vi;

end