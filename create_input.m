function[tis, Cas, Lis, vis, ts] = create_input(tiso, dTt, dTc, ISI, Ca, N)

if nargin < 6
    N = 500;
end

%% create idealized input
Ts = cumsum([tiso-(dTt*3+dTc*2+ISI); dTc; dTc; ISI; dTt; dTt; dTt]);
tx = linspace(0, Ts(end), 10*N);

% specify velocity
vx = nan(size(tx));
vx(tx <= Ts(1)) = 0;
vx(tx > Ts(1) & tx <=  Ts(2)) = .4545;
vx(tx > Ts(2) & tx <=  Ts(3)) = -.4545;
vx(tx > Ts(3) & tx <=  Ts(4)) = 0;
vx(tx > Ts(4) & tx <=  Ts(5)) = .4545;
vx(tx > Ts(5) & tx <=  Ts(6)) = 0;
vx(tx > Ts(6) & tx <=  Ts(7)) = -.4545;

% integrate to get length
Lx = cumtrapz(tx, vx);

% interpolate
ti = linspace(0, Ts(end), N);
vi = interp1(tx, vx, ti);
Li = interp1(tx, Lx, ti);

%% combine different pCa trials
vis = [];
Cas = [];
Lis = [];

for k = 1:length(Ca)
    vis = [vis vi];
    Cas = [Cas Ca(k) * ones(size(vi))];
    Lis = [Lis Li];
end

tis = linspace(0, Ts(end)*length(Ca), length(ti)*length(Ca));
ts = 0:tiso:(tiso*(length(Ca)-1));

end

