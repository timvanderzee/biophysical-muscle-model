function[IG] = get_initial_guess(Xdata, parms, visualize)

if ~isfield(Xdata, 'idA')
    Xdata.idA = 1:length(Xdata.Cas);
end

parms.ti = Xdata.t;
parms.vts = Xdata.v;
parms.Cas = Xdata.Cas;
parms.Lts = Xdata.L * parms.gamma;

tis = parms.ti(Xdata.idA);
Cas = parms.Cas(Xdata.idA);
vts = parms.vts(Xdata.idA);
Lts = parms.Lts(Xdata.idA);

nparms = parms;
nparms.vts = [0 0];
nparms.ti = [0 5];
nparms.Lts = [0 0];
nparms.Cas = Cas(1) * [1 1];

nparms.PE_isw_SE = 1;

% find L0
Lcost = @(dlse, parms) ((parms.kse0*(exp(parms.kse*dlse)-1)) - (parms.kpe * (0 - dlse - parms.Lce0))).^2;

L0 = fminsearch(@(L) Lcost(L, nparms), 10);
Fse0 = nparms.Fse_func(L0, nparms);


% [Q00, Q20, lce0, Q10, Fse0] = find_steady_state(Q0, p0, q0, parms, 'regular');
x0 = [0 0 Fse0 0 1 0]';

odeopt = odeset('maxstep', 1e-2);

% isometric contraction
sol0 = ode15s(@(t,y) fiber_dynamics_explicit_length_v2(t,y, nparms), [0 5], x0, odeopt);

% dynamic contraction
nparms.vts = vts;
nparms.Lts = Lts;
nparms.ti = tis;
nparms.Cas = Cas;

sol = ode15s(@(t,y) fiber_dynamics_explicit_length_v2(t,y,nparms), [0 max(tis)], sol0.y(:,end), odeopt);

%% extract the output and interpolate
% states
Q0i     = interp1(sol.x, sol.y(1,:), tis); % zero-order moment
Q2i     = interp1(sol.x, sol.y(2,:), tis); % second-order moment
Fsei    = interp1(sol.x, sol.y(3,:), tis); % length
Noni    = interp1(sol.x, sol.y(4,:), tis); % thin filament activation
DRXi    = interp1(sol.x, sol.y(5,:), tis); % thick filament activation
Ri      = interp1(sol.x, sol.y(6,:), tis); %

%% determine Fpe
dLse = max(nparms.Lse_func(Fsei, nparms), 0); % can't be negative
Lce = Lts - dLse;

if isfield(nparms, 'PE_isw_SE') % PE in series with SE
    dLce = Lce - nparms.Lce0;

    if nparms.K*dLce < 10
        kpe = nparms.kpe .* (1 - 1./(exp(nparms.K*dLce)+1));
        Fpe = nparms.kpe * log(1+exp(dLce*nparms.K))/nparms.K;
    else
        kpe = nparms.kpe;
        Fpe = nparms.kpe * dLce;
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
k1 = [nparms.k11 nparms.k12];
k2 = [nparms.k21 -nparms.k22];

% get functions
[IGf, IGef] = get_IG_IGEf(nparms.approx);

[Q0dot, Q1dot] = CrossBridge_Dynamics(Q0i, pi, qi, nparms.f, nparms.w, k1, k2, IGef, Noni, DRXi, IGf, nparms.b, nparms.k, Ri, nparms.dLcrit, nparms.ps2);

% velocity - independent derivative
F0dot  = Q1dot + Q0dot;

kse = nparms.kse * (Fsei + nparms.kse0);
Ldi = (vts .* nparms.gamma .* kse - F0dot) ./  (Q0i + kse + kpe);
IG.Ldi = Ldi;

% other variables
IG.Q1i = Q1i;
IG.Fcei = Fce;
IG.Fpei = Fpe;

IG.pi = pi;
IG.qi = qi;
IG.Q00i = Q00i;

%% passive
Lcost = @(dlse, Lis, parms) ((parms.kse0*(exp(parms.kse*dlse)-1)) - (parms.kpe * (Lis - dlse - parms.Lce0))).^2;

for i = 1:length(parms.Lts)
    dlse(i) = fminsearch(@(L) Lcost(L, parms.Lts(i), nparms), 0);
end

dLce = (parms.Lts - dlse) - nparms.Lce0;

Fse = nparms.Fse_func(dlse, nparms) * nparms.Fscale;
Fpe = nparms.kpe * dLce  * nparms.Fscale;

if visualize
figure(1)
subplot(414)
plot(Xdata.t, Fse, '-', Xdata.t, Fpe, '--')
end

% initial guess
IG.F9 = Fse;

%% plotting
oFi = (IG.Fsei) * nparms.Fscale;

if visualize
%     figure(1r
    subplot(414); hold on
    
    plot(Xdata.t(Xdata.idA), oFi,'b'); hold on

    
    if isfield(Xdata, 'idF') 
        plot(Xdata.t(Xdata.idF), oFi(Xdata.idF),'m*', 'markersize', 1); hold on
    end
    
    if isfield(Xdata, 'idC')
        plot(Xdata.t(Xdata.idC), oFi(Xdata.idC),'g.'); hold on
    end
end

end
