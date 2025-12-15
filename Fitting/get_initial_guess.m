function[IG] = get_initial_guess(tis, Cas, vts, Lts, parms)

parms.vts = [0 0];
parms.ti = [0 5];
parms.Lts = [0 0];
parms.Cas = Cas(1) * [1 1];

parms.PE_isw_SE = 1;

% find L0
Lcost = @(dlse, parms) ((parms.kse0*(exp(parms.kse*dlse)-1)) - (parms.kpe * (0 - dlse - parms.Lce0))).^2;

L0 = fminsearch(@(L) Lcost(L, parms), 10);
Fse0 = parms.Fse_func(L0, parms);


% [Q00, Q20, lce0, Q10, Fse0] = find_steady_state(Q0, p0, q0, parms, 'regular');
x0 = [0 0 Fse0 0 1 0]';

odeopt = odeset('maxstep', 1e-2);

% isometric contraction
sol0 = ode15s(@(t,y) fiber_dynamics_explicit_length_v2(t,y, parms), [0 5], x0, odeopt);

% dynamic contraction
parms.vts = vts;
parms.Lts = Lts;
parms.ti = tis;
parms.Cas = Cas;

sol = ode15s(@(t,y) fiber_dynamics_explicit_length_v2(t,y,parms), [0 max(tis)], sol0.y(:,end), odeopt);

%% extract the output and interpolate
% states
Q0i     = interp1(sol.x, sol.y(1,:), tis); % zero-order moment
Q2i     = interp1(sol.x, sol.y(2,:), tis); % second-order moment
Fsei    = interp1(sol.x, sol.y(3,:), tis); % length
Noni    = interp1(sol.x, sol.y(4,:), tis); % thin filament activation
DRXi    = interp1(sol.x, sol.y(5,:), tis); % thick filament activation
Ri      = interp1(sol.x, sol.y(6,:), tis); %

%% determine Fpe
dLse = max(parms.Lse_func(Fsei, parms), 0); % can't be negative
Lce = Lts - dLse;

if isfield(parms, 'PE_isw_SE') % PE in series with SE
    dLce = Lce - parms.Lce0;

    if parms.K*dLce < 10
        kpe = parms.kpe .* (1 - 1./(exp(parms.K*dLce)+1));
        Fpe = parms.kpe * log(1+exp(dLce*parms.K))/parms.K;
    else
        kpe = parms.kpe;
        Fpe = parms.kpe * dLce;
    end
end

%% Determine cross-bridge states
Fce = Fsei - Fpe;
Q1i = Fce - Q0i;

% mean and standard deviation
Q00i = max(Q0i, 1e-3);
pi = Q1i./Q00i;
qi = Q2i./Q00i - pi.^2;
qi = max(qi, 0);

% figure(1)
% plot(qi)

%% save to struct
% states
IG.Q0i = Q0i;
IG.Q2i = Q2i;
IG.Fsei = Fsei;
IG.Noni = Noni;
IG.DRXi = DRXi;
IG.Ri = Ri;
%     IG.Lti = Lti;

%% compute Ldi
% kpe = 0;

% points where integrals is evaluated
k1 = [parms.k11 parms.k12];
k2 = [parms.k21 -parms.k22];

% get functions
[IGf, IGef] = get_IG_IGEf(parms.approx);

[Q0dot, Q1dot] = CrossBridge_Dynamics(Q0i, pi, qi, parms.f, parms.w, k1, k2, IGef, Noni, DRXi, IGf, parms.b, parms.k, Ri, parms.dLcrit, parms.ps2);

% velocity - independent derivative
F0dot  = Q1dot + Q0dot;

kse = parms.kse * (Fsei + parms.kse0);
Ldi = (vts .* parms.gamma .* kse - F0dot) ./  (Q0i + kse + kpe);
IG.Ldi = Ldi;

% other variables
IG.Q1i = Q1i;
IG.Fcei = Fce;
IG.Fpei = Fpe;

IG.pi = pi;
IG.qi = qi;
IG.Q00i = Q00i;

end
