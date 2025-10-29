function[yp] = fiber_dynamics_explicit_no_tendon_full(t,y, parms)

% Get velocity and calcium
if numel(parms.vts) == 1
    vMtilda = parms.vts;
%     lMtilda = parms.Lts;
else
    vMtilda = interp1(parms.ti, parms.vts, t);
%     lMtilda = interp1(parms.ti, parms.Lts, t);
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

R = 0;
dRdt = 0;

% displacement from start
xi = parms.xi + (L - parms.lce0);

% compute stiffness
n(n<0) = 0;
Q0 = trapz(xi(:), n);
Q1 = trapz(xi(:), n.*xi(:));
F = Q0 + Q1;

% avoid small numbers
k   = parms.K;
F   = log(1+exp(F*k))/k;
% Q00   = log(1+exp(Q0*k))/k;
% Q00 = Q0;

% points where integrals is evaluated
k1 = [parms.k11 parms.k12];
k2 = [parms.k21 -parms.k22];

% attachment and detachment at each strain
beta = parms.f_func(xi, parms.f, parms.w);
phi = -(parms.g_func(xi, k1(1), -k1(2)) + parms.g_func(xi, k2(1), -k2(2))) .* n';   

% change in cross-bridge attachment
ndot = DRX * (beta * (Non - Q0)) + phi;

% first determine contraction velocity
Qdot = trapz(xi(:), [ndot(:) xi(:).*ndot(:)]);
F0dot  = Qdot(1) + Qdot(2);

[Jon, Joff] = ThinFilament_Dynamics(Act, Q0, Non, parms.kon, parms.koff, parms.koop, 1);
dNondt = Jon - Joff;

[J1, J2] = ThickFilament_Dynamics(Q0, F, DRX, parms.J1, parms.J2, parms.JF, 1, R);

% note: the term "Q0dot + Rdot" is the part of Q0dot due to regular
% attachment. explanation: if Rdot is great, it means that Q0 is losing to
% R. this implies that for a given Q0dot, the part due to regular
% attachment must be greater. thus, we need to add Rdot to Q0dot
dDRXdt = J1 - J2 - (Qdot(1) + dRdt);

kse = parms.kse * (F + parms.kse0);

Ld  = (vMtilda .* parms.gamma .* kse - F0dot) ./ (Q0 + kse);

yp = [ndot(:); Ld; dNondt; dDRXdt];

end