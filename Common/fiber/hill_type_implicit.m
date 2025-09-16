function[error, Fse] = hill_type_implicit(t, x, xdot, parms)

gamma = (.5*parms.s)/parms.h;

lmtc = interp1(parms.ti, parms.Lts, t) * gamma;

if numel(parms.Cas, 1)
    Ca = parms.Cas;
else
    Ca      = interp1(parms.ti, parms.Cas, t);
end
% Ca   = interp1(parms.ti, parms.Cas, t);

% states
lce = x(1);
%     lmtc = x(2);

vce = xdot(1);
%     vmtc = xdot(2);

Fce_rel = parms.vF_func(vce, parms);

% activation from Ca
a = parms.actfunc(Ca, parms);
a(a<parms.amin) = parms.amin;

Fce = a * Fce_rel;

% elastic elements
Lse = lmtc - lce;

if isfield(parms, 'Fpece_func')
    Fpe = parms.Fpece_func(lce, parms);
else
    Fpe = 0;
end

Fse = parms.Fse_func(Lse, parms);
Fse(Lse < 0) = 0;

% error terms
error(1,1) = Fse - Fce - Fpe;
%     error(2,1) = parms.vmtc - vmtc;

end