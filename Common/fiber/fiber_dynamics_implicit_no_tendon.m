function[error] = fiber_dynamics_implicit_no_tendon(t,y,yp, parms)

% Get velocity and calcium
vMtilda = interp1(parms.ti, parms.vts, t);

if numel(parms.Cas, 1)
    Ca = parms.Cas;
else
    Ca      = interp1(parms.ti, parms.Cas, t);
end

% States
Q0  = y(1);
Q1  = y(2);
Q2  = y(3);
Non = y(4);
DRX = y(5);
Ld  = y(6);

Act = parms.actfunc(Ca, parms);

% in case we don't do thin filament cooperative activation
if (parms.kon == 0) && (parms.koff == 0) && (parms.koop == 0)
    Non = Act;
end

% State derivatives
dQ0dt  = yp(1);
dQ1dt  = yp(2);
dQ2dt  = yp(3);
dNondt = yp(4);
dDRXdt = yp(5);

% ripped
if length(y) > 6
    R   = y(7);
    dRdt = yp(7);
else
    R = 0;
    dRdt = 0;
end

% Cross-bridge states
k = 100;
F = Q1 + Q0;
F   = log(1+exp(F*k))/k;
Q00 = log(1+exp(Q0*k))/k;

% mean and standard deviation
p = Q1./Q00; 
q = Q2./Q00 - p.^2;  
q = log(1+exp(q*k))/k;

% Thin and thick filament
[error_thin, ~] = ThinEquilibrium(Act, Q0, Non, dNondt, parms.kon, parms.koff, parms.koop, parms.act * parms.Noverlap);
[error_thick, ~] = ThickEquilibrium(Q0, dQ0dt, F, DRX, dDRXdt, parms.J1, parms.J2, parms.JF, parms.act * parms.Noverlap, R);

% Cross-bridge dynamics
[error_Q0, error_Q1, error_Q2, error_R, Fdot] = MuscleEquilibrium(Q0, Q1, p, q, dQ0dt, dQ1dt, dQ2dt, parms.f, parms.w, parms.k11, parms.k12, parms.k21, parms.k22, Non, Ld, DRX, dRdt, parms.b, parms.k, R, parms.dLcrit);

% Length dynamics
[error_length] = LengthEquilibrium(Q0, F, Fdot, Ld, vMtilda, parms.kse0, parms.kse);

% Combined error
if length(y) > 6
error = [error_thin; error_thick; error_Q0; error_Q1; error_Q2; error_R; error_length];
else
error = [error_thin; error_thick; error_Q0; error_Q1; error_Q2; error_length];
end

end