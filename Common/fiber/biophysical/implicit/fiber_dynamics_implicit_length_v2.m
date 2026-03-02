function[error] = fiber_dynamics_implicit_length_v2(t,y,yp, parms)

% in this alternative version we have Fse as a state, which may be less prone
% to error

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
Q2  = y(2);
Fse  = y(3);
Non = y(4);
DRX = y(5);
R = y(6);

if length(y) > 6
    Lts = y(7);
    dLdt = yp(7);
else
    if numel(parms.Lts) == 1
        Lts = parms.Lts;
    else
        Lts = interp1(parms.ti, parms.Lts, t);
    end
end

% SE
Fse = max(Fse, 1e-6);
dLse = parms.Lse_func(Fse, parms);
kse = parms.kse * (Fse + parms.kse0);

% PE
Lce = Lts - dLse;
kpe = 0;
Fpe = 0;

if isfield(parms, 'PE_isw_SE') % PE in series with SE
    if parms.PE_isw_SE
        kpe = parms.kpe_func(Lce, parms);
        Fpe = parms.Fpe_func(Lce, parms);
    end
end

% CE
Fce = Fse - Fpe;

% States
dQ0dt  = yp(1);
dQ2dt  = yp(2);
dFsedt  = yp(3);
dNondt = yp(4);
dDRXdt = yp(5);
dRdt = yp(6);

% mean and standard deviation
% k   = parms.K;
% Q00 = log(1+exp(Q0*k))/k; % note: goes to inf for large k, may need another function
Q00 = max(Q0, 1e-6);
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

% points where integrals is evaluated
k1 = [parms.k11 parms.k12];
k2 = [parms.k21 -parms.k22];

% get functions
[IG, IGef] = get_IG_IGEf(parms.approx);

% Compute Qdot
[Q0dot, Q1dot, Q2dot, Rdot] = CrossBridge_Dynamics(Q0, p, q, parms.f, parms.w, k1, k2, IGef, Non, DRX, IG, parms.b, parms.k, R, parms.dLcrit, parms.ps2);

% velocity - independent derivative
F0dot  = Q1dot + Q0dot;

% Calculate Ldot
Ld = (vMtilda .* parms.gamma .* kse - F0dot) /  (Q0 + kse + kpe);

% Cross-bridge errors
error_Q0 = dQ0dt - Q0dot;
error_Q2 = dQ2dt - (Q2dot + 2 * Ld .* Q1);

% Ripped dynamics
error_R = dRdt - Rdot;

% force error?
% error_force = Fse - Fce - Fpe;

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