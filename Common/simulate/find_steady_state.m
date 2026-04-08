function[Q00, Q20, lce0, Q10, Fse, Fpe, Fce] = find_steady_state(Q0, p0, q0, L0, parms, type)
% finds the steady-state length by solving Fce + Fpe = Fse
% also returns Q00 (obeying log function) and Q20 (from q0, p0 and Q0)

% avoid very small numbers
% Q00 = log(1+exp(Q0))/parms.K; % note: goes to inf for large k, may need another function
Q00 = max(Q0, 1e-6);

% provide Q2
Q20 = (q0 + p0^2) * Q00;

% compute Fce
Q10 = Q00 * p0;
Fce = (Q10 + Q00);


% compute the length that satisfies force constraint
if strcmp(type, 'regular')
    costfunc = @(L, Fce, L0, parms) (Fce + (parms.kpe*(L-parms.Lce0).*(L>parms.Lce0)) - (parms.kse0*(exp(parms.kse*(L0-L))-1))) .^2;
elseif strcmp(type, 'adjusted')
    costfunc = @(L, Fce, L0, parms) (Fce - (parms.kse0*(exp(parms.kse*-L)-1))) .^2;
end

lce0 = fminsearch(@(L) costfunc(L, Fce, L0, parms), 0);

if strcmp(type, 'regular')
    Fpe = parms.kpe*(lce0-parms.Lce0).*(lce0>parms.Lce0);
else
    Fpe = 0;
end

Fse0 = parms.kse0*(exp(parms.kse*(L0-lce0))-1);

Fse = Fce + Fpe;

end