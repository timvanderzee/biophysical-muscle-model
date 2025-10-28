function[error] = fiber_dynamics_implicit_no_tendon_full(t,y,yp, parms)

% Get velocity and calcium
if numel(parms.vts) == 1
    vMtilda = parms.vts;
    lMtilda = parms.Lts;
else
    vMtilda = interp1(parms.ti, parms.vts, t);
    lMtilda = interp1(parms.ti, parms.Lts, t);
end

if numel(parms.Cas) == 1
    Ca = parms.Cas;
else
    Ca = interp1(parms.ti, parms.Cas, t);
end

% activation
Act = parms.actfunc(Ca, parms);

% States
n   = y(1:length(parms.xi));
L   = y(end-2);
Non = y(end-1);
DRX = y(end-0);

% State derivatives
dndt  = yp(1:length(parms.xi));
Ld     = yp(end-2);
dNondt = yp(end-1);
dDRXdt = yp(end-0);

R = 0;
dRdt = 0;

% calculate force
dlse = lMtilda - L;
F = parms.Fse_func(dlse, parms);

% displacement from start
xi = parms.xi + (L - parms.lce0);

% compute stiffness
n(n<0) = 0;
Q0 = trapz(xi(:), n);

% Cross-bridge dynamics
[error_n, Qdot] = MuscleEquilibrium_full(n, dndt, Q0, Non, DRX, parms.f, parms.w, xi, parms.k11, parms.k12, parms.k21, parms.k22, parms.f_func, parms.g_func);
dQ0dt = Qdot(1);
F0dot  = dQ0dt + Qdot(2);

% Thin and thick filament
if (parms.kon == 0) && (parms.koff == 0) && (parms.koop == 0)
    error_thin = dNondt - ((Act - Non) / .005);
else 
    [error_thin, ~] = ThinEquilibrium(Act, Q0, Non, dNondt, parms.kon, parms.koff, parms.koop, parms.act * parms.Noverlap);
end

[error_thick, ~] = ThickEquilibrium(Q0, dQ0dt, F, DRX, dDRXdt, parms.J1, parms.J2, parms.JF, parms.act * parms.Noverlap, R, dRdt);

% Length dynamics
[error_length] = LengthEquilibrium(Q0, F, F0dot, Ld, vMtilda, parms.kse0, parms.kse, parms.gamma);

% Combined error
if length(y) < 6
    error_thick = [];
end

error = [error_n; error_length; error_thin; error_thick];

end