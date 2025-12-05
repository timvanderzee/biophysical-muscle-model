function[error] = fiber_dynamics_implicit_length(t,y,yp, parms)

% Get velocity and calcium
if numel(parms.vts) == 1
    vMtilda = parms.vts;
else
    vMtilda = interp1(parms.ti, parms.vts, t);
end

% if numel(parms.Lts) == 1
%     Lts = parms.Lts;
% else
%     Lts = interp1(parms.ti, parms.Lts, t);
% end

if numel(parms.Cas) == 1
    Ca = parms.Cas;
else
    Ca = interp1(parms.ti, parms.Cas, t);
end

% activation
Act = parms.actfunc(Ca, parms);

% States
Q0  = y(1);
Q2  = y(2);
Lce  = y(3);
Non = y(4);
DRX = y(5);
R = y(6);
Lts = y(7);

% PE
kpe = parms.kpe_func(Lce, parms);
Fpe = parms.Fpe_func(Lce, parms);

% SE
dLse = Lts - Lce;
Fse = parms.Fse_func(dLse, parms);
kse = parms.kse_func(dLse, parms);
% kse = parms.kse * (Fse .* (Fse > 0) + parms.kse0);

% CE
Fce = Fse - Fpe;

% States
dQ0dt  = yp(1);
dQ2dt  = yp(2);
dLcedt  = yp(3);
dNondt = yp(4);
dDRXdt = yp(5);
dRdt = yp(6);
dLdt = yp(7);

% mean and standard deviation
% k   = parms.K;
% Q00 = log(1+exp(Q0*k))/k; % note: goes to inf for large k, may need another function
Q00 = max(Q0, 1e-4);
Q1 = Fce - Q00;
p = Q1./Q00; 
q = Q2./Q00 - p.^2;  
% q = log(1+exp(q*k))/k;
q = max(q, 1e-6);

% Thin and thick filament
if (parms.kon == 0) && (parms.koff == 0) && (parms.koop == 0)
    error_thin = dNondt - ((Act - Non) / .005);
else 
    [error_thin, ~] = ThinEquilibrium(Act, Q0, Non, dNondt, parms.kon, parms.koff, parms.koop, parms.act * parms.Noverlap);
end

[error_thick, ~] = ThickEquilibrium(Q0, dQ0dt, Fce, DRX, dDRXdt, parms.J1, parms.J2, parms.JF, parms.act * parms.Noverlap, R, dRdt);

% Cross-bridge dynamics
[error_Q0, ~, error_Q2, F0dot, Rdot] = MuscleEquilibrium(Q0, Q1, p, q, dQ0dt, 0, dQ2dt, parms.f, parms.w, parms.k11, parms.k12, parms.k21, parms.k22, Non, dLcedt, DRX, parms.b, parms.k, R, parms.dLcrit, parms.ps2, parms.approx);

% Rippid dynamics
error_R = dRdt - Rdot;

% Length dynamics
% error_length = dLcedt  .* (Q0 + kse + kpe) - (vMtilda .* parms.gamma .* kse - F0dot);
error_length = LengthEquilibrium(Q0, F0dot, dLcedt, vMtilda, kse, kpe, parms.gamma);

error_vel = dLdt - vMtilda*parms.gamma;

% Combined error
error = [error_Q0; error_Q2; error_length; error_thin; error_thick; error_R; error_vel];

end