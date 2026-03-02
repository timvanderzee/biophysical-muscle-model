function[yp, Fce, Q0] = fiber_dynamics_explicit_discretized(t,y, parms)

% get input
[vMtilda, ~, Ca] = get_input_from_parms(t, parms);

% States
n   = y(1:length(parms.xi));
L   = y(end-3);
Non = y(end-2);
DRX = y(end-1);
R   = y(end-0);

% displacement from start
xi = parms.xi + (L - parms.lce0);

% overwrite if not using MoC
if isfield(parms, 'method_of_characteristics')
    if ~parms.method_of_characteristics
        xi = parms.xi;
    end
end

% compute stiffness
n(n<0) = 0;
Q0 = trapz(xi(:), n);

% if isfield(parms, 'method_of_characteristics')
%     if parms.method_of_characteristics
        Q1 = trapz(xi(:), n.*xi(:));
        Fce = Q0 + Q1;
%     end

% % overwrite if not using MoC
% if isfield(parms, 'method_of_characteristics')
%     if ~parms.method_of_characteristics
%         F = parms.F;
%     end
% end

% avoid small numbers
k   = parms.K;
Fce   = log(1+exp(Fce*k))/k;

% points where integrals is evaluated
k1 = [parms.k11 parms.k12];
k2 = [parms.k21 -parms.k22];

% attachment and detachment at each strain
beta = parms.f_func(xi, parms.f, parms.w, 0);
phi = -(parms.g_func(xi, k1(1), -k1(2)) + parms.g_func(xi, k2(1), -k2(2))) .* n';   

% forcible detachment
gamma = parms.f_func(xi, parms.b, parms.w, 0) * R;
% phiR = -parms.k*(xi > (parms.dLcrit)) .* n' + gamma;
phiR = -parms.k * (1/2 * (tanh((xi - parms.dLcrit)*20) + 1)) .* n' + gamma;
dRdt = -trapz(xi, phiR);

% total phi
phiT = phi + phiR;

% change in cross-bridge attachment
ndot = DRX * (beta * (Non - Q0)) + phiT;

% first determine contraction velocity
Qdot = trapz(xi(:), [ndot(:) xi(:).*ndot(:)]);
F0dot  = Qdot(1) + Qdot(2);

if (parms.kon == 0) && (parms.koff == 0) && (parms.koop == 0)
	Act = Ca.^parms.n./(parms.kappa^parms.n+Ca.^parms.n); % sigmoidal function
    dNondt = ((Act - Non) / .005);
else
    [Jon, Joff] = ThinFilament_Dynamics(Ca, Q0, Non, parms.kon, parms.koff, parms.koop, 1);
    dNondt = Jon - Joff;
end

[J1, J2] = ThickFilament_Dynamics(Q0, Fce, DRX, parms.J1, parms.J2, parms.JF, 1, R);

% note: the term "Q0dot + Rdot" is the part of Q0dot due to regular
% attachment. explanation: if Rdot is great, it means that Q0 is losing to
% R. this implies that for a given Q0dot, the part due to regular
% attachment must be greater. thus, we need to add Rdot to Q0dot
dDRXdt = J1 - J2 - (Qdot(1) + dRdt);

% compute velocity
[Fpe, kpe] = get_PE(L, parms);
Fse = Fce + Fpe;
kse = parms.kse * (Fse + parms.kse0);
Ld  = (vMtilda .* parms.gamma .* kse - F0dot) ./ (Q0 + kse + kpe);

% total state derivative
yp = [ndot(:); Ld; dNondt; dDRXdt; dRdt];

end