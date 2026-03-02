function[error] = fiber_dynamics_implicit_approximated(t,y,yp, parms)

% get input
[vMtilda, Lts, Ca] = get_input_from_parms(t, parms);
    
% States
Q0      = y(1);
Q2      = y(2);
Fse     = y(3);
Non     = y(4);
DRX     = y(5);
R       = y(6);

% State derivatives
dQ0dt   = yp(1);
dQ2dt   = yp(2);
dFsedt  = yp(3);
dNondt  = yp(4);
dDRXdt  = yp(5);
dRdt    = yp(6);

% get PE and SE stiffness, CE force
[Fce, kpe, kse] = get_force_and_stiffness(Fse, Lts, parms);

% get distribution mean and standard deviation
[p, q, Q1] = get_pq(Q0, Q2, Fce);

% Thin and thick filament
if (parms.kon == 0) && (parms.koff == 0) && (parms.koop == 0)
    Act = Ca.^parms.n./(parms.kappa^parms.n+Ca.^parms.n); % sigmoidal function
    error_thin = dNondt - ((Act - Non) / .005);
else 
    [error_thin, ~] = ThinEquilibrium(Ca, Q0, Non, dNondt, parms.kon, parms.koff, parms.koop, parms.act * parms.Noverlap);
end

[error_thick, ~] = ThickEquilibrium(Q0, dQ0dt, Fce, DRX, dDRXdt, parms.J1, parms.J2, parms.JF, parms.act * parms.Noverlap, R, dRdt);

% points where integrals is evaluated
k1 = [parms.k11 parms.k12];
k2 = [parms.k21 -parms.k22];

% get functions
[IG, IGef] = get_IG_IGEf(parms.approx);

% Compute Qdot
[Q0dot, Q1dot, Q2dot, Rdot] = CrossBridge_Dynamics(Q0, p, q, parms.f, parms.w, k1, k2, IGef, Non, DRX, IG, parms.b, parms.k, R, parms.dLcrit, 0);

% velocity - independent derivative
F0dot  = Q1dot + Q0dot;

% Calculate Ldot
Ld = (vMtilda .* parms.gamma .* kse - F0dot) /  (Q0 + kse + kpe);

% Cross-bridge errors
error_Q0 = dQ0dt - Q0dot;
error_Q2 = dQ2dt - (Q2dot + 2 * Ld .* Q1);

% Ripped dynamics
error_R = dRdt - Rdot;

% force-rate error
dQ1dt = Q1dot + 1 * Ld .* Q0;
dFcedt = dQ1dt + dQ0dt;
dFpedt = Ld * kpe;
error_force = dFsedt - dFcedt - dFpedt;

if length(y) > 6
    error_vel = dLdt - vMtilda*parms.gamma;
else
    error_vel = [];
end

% Combined error
error = [error_Q0; error_Q2; error_force; error_thin; error_thick; error_R; error_vel];

end