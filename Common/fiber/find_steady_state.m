function[Q00, Q20, lce0] = find_steady_state(Q0, p0, q0, parms)
% finds the steady-state length by solving Fce + Fpe = Fse
% also returns Q00 (obeying log function) and Q20 (from q0, p0 and Q0)

% avoid very small numbers
Q00 = log(1+exp(Q0))/parms.K; % note: goes to inf for large k, may need another function

% provide Q2
Q20 = (q0 + p0^2) * Q00;

% compute Fce
Q10 = Q00 * p0;
Fce = (Q10 + Q00);

% compute the length that satisfies force constraint
costfunc = @(L, Fce, parms) (Fce + (parms.kpe*(L-parms.lmtc0).*(L>parms.lmtc0)+parms.Fpe0) - (parms.kse0*(exp(parms.kse*-L)-1))) .^2;
lce0 = fminsearch(@(L) costfunc(L, Fce, parms), 0);

end