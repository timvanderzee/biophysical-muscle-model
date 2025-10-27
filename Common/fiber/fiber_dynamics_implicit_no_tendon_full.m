function[error] = fiber_dynamics_implicit_no_tendon_full(t,y,yp, parms)

% Get velocity and calcium
if numel(parms.vts) == 1
    vMtilda = parms.vts;
else
    vMtilda = interp1(parms.ti, parms.vts, t);
end

if numel(parms.Cas) == 1
    Ca = parms.Cas;
else
    Ca = interp1(parms.ti, parms.Cas, t);
end

% activation
Act = parms.actfunc(Ca, parms);

% States
Q0  = y(1);
Q1  = y(2);
Q2  = y(3);
Ld  = y(4);
Non = y(5);

% State derivatives
dQ0dt  = yp(1);
dQ1dt  = yp(2);
dQ2dt  = yp(3);

dNondt = yp(5);

if length(y) > 5
%     Non = y(5);
    DRX = y(6);
%     dNondt = yp(5);
    dDRXdt = yp(6);
else
%     Non = Act;
    DRX = 1 - Q0;
%     dNondt = 0;
    dDRXdt = 0;
end

% % in case we don't do thin filament cooperative activation
% if (parms.kon == 0) && (parms.koff == 0) && (parms.koop == 0)
%     Non = Act;
% end

% ripped
if length(y) > 6
    R   = y(7);
    dRdt = yp(7);
else
    R = 0;
    dRdt = 0;
end

% Cross-bridge states
k   = parms.K;
F   = Q1 + Q0;
F   = log(1+exp(F*k))/k;
Q00 = log(1+exp(Q0*k))/k;

% mean and standard deviation
p = Q1./Q00; 
q = Q2./Q00 - p.^2;  
q = log(1+exp(q*k))/k;

% Thin and thick filament
if (parms.kon == 0) && (parms.koff == 0) && (parms.koop == 0)
    error_thin = dNondt - ((Act - Non) / .005);
else 
    [error_thin, ~] = ThinEquilibrium(Act, Q0, Non, dNondt, parms.kon, parms.koff, parms.koop, parms.act * parms.Noverlap);
end

[error_thick, ~] = ThickEquilibrium(Q0, dQ0dt, F, DRX, dDRXdt, parms.J1, parms.J2, parms.JF, parms.act * parms.Noverlap, R, dRdt);

% Cross-bridge dynamics
% [error_Q0, error_Q1, error_Q2, error_R, F0dot] = MuscleEquilibrium(Q0, Q1, p, q, dQ0dt, dQ1dt, dQ2dt, parms.f, parms.w, parms.k11, parms.k12, parms.k21, parms.k22, Non, Ld, DRX, dRdt, parms.b, parms.k, R, parms.dLcrit, parms.ps2, parms.approx);
[error_x, error_R, F0dot] = MuscleEquilibrium_full(x, dQ0dt, dQ1dt, dQ2dt, parms.f, parms.w, parms.k11, parms.k12, parms.k21, parms.k22, Non, Ld, DRX, dRdt, parms.b, parms.k, R, parms.dLcrit, parms.ps2, parms.approx);

% Length dynamics
[error_length] = LengthEquilibrium(Q0, F, F0dot, Ld, vMtilda, parms.kse0, parms.kse, parms.gamma);

% Combined error
if length(y) < 6
%     error_thin = [];
    error_thick = [];
end

if length(y) < 7
    error_R = [];
end

error = [error_Q0; error_Q1; error_Q2; error_length; error_thin; error_thick; error_R];

end