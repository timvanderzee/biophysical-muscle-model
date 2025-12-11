function[dx, Fce, Q0] = fiber_dynamics_explicit_length_v2(t,y, parms)

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
Fse  = y(3);
Non = y(4);
DRX = y(5);
R = y(6);

if length(y) > 6
    Lts = y(7);
%     dLdt = yp(7);
else
    if numel(parms.Lts) == 1
        Lts = parms.Lts;
    else
        Lts = interp1(parms.ti, parms.Lts, t);
    end
end

% SE
dLse = max(parms.Lse_func(Fse, parms), 0); % can't be negative
kse = parms.kse * (Fse + parms.kse0);

% PE
Lce = Lts - dLse;
kpe = 0;
Fpe = 0;
% 
% if isfield(parms, 'PE_isw_SE') % PE in series with SE
%     if parms.PE_isw_SE
%         kpe = parms.kpe_func(Lce, parms);
%         Fpe = parms.Fpe_func(Lce, parms);
%     end
% end
if isfield(parms, 'PE_isw_SE') % PE in series with SE
    dLce = Lce - parms.Lce0;

    if parms.K*dLce < 10
%         kpe = parms.kpe .* (1 - 1./(exp(parms.K*dLce)+1));
        Fpe = parms.kpe * log(1+exp(dLce*parms.K))/parms.K;
    else
%         kpe = parms.kpe;
        Fpe = parms.kpe * dLce;
    end
end

% CE
Fce = Fse - Fpe;
Fce = max(Fce, 0);

% mean and standard deviation
% k   = parms.K;
% Q00 = log(1+exp(Q0*k))/k; % note: goes to inf for large k, may need another function
Q00 = max(Q0, 1e-6);
Q1 = Fce - Q00;
p = Q1./Q00; 
q = Q2./Q00 - p.^2;  
% q = log(1+exp(q*k))/k;
q = max(q, 1e-6);

if isfield(parms, 'FL_overlap')
    if parms.FL_overlap
        h = (.5 * parms.s / parms.gamma); % powerstroke size
        L_hs = Lce * h * 1e9 + 1.3e3; % [nm]
        Ntot = return_f_overlap(L_hs, parms);
    else
        Ntot = 1;
    end
else
    Ntot = 1;
end

% Thin and thick filament
if (parms.kon == 0) && (parms.koff == 0) && (parms.koop == 0)
    dNondt = Ntot * ((Act - Non) / .005);
else 
    [Jon, Joff] = ThinFilament_Dynamics(Act, Q0, Non, parms.kon, parms.koff, parms.koop, Ntot);
    dNondt = Jon - Joff;
end

[J1, J2] = ThickFilament_Dynamics(Q0, Fce, DRX, parms.J1, parms.J2, parms.JF, 1, R);

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

% shift the cross-bridge distribution
dQ0dt = Q0dot;
dQ2dt = (Q2dot + 2 * Ld .* Q1);

% force-rate error
dQ1dt = Q1dot + 1 * Ld .* Q0;
dFcedt = dQ1dt + dQ0dt;
dFpedt = Ld * kpe;
dFsedt = dFcedt + dFpedt;

% note: the term "Q0dot + Rdot" is the part of Q0dot due to regular
% attachment. explanation: if Rdot is great, it means that Q0 is losing to
% R. this implies that for a given Q0dot, the part due to regular
% attachment must be greater. thus, we need to add Rdot to Q0dot
dDRXdt = J1 - J2 - (Q0dot + Rdot);

if length(y) > 6
    Ldot = vMtilda .* parms.gamma;
else
    Ldot = [];
end

dx = [dQ0dt; dQ2dt; dFsedt; dNondt; dDRXdt; Rdot; Ldot];

end