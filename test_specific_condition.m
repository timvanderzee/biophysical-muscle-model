clear all; close all; clc

iFs = [1 2 3, 5, 6, 7, 8, 10, 11];
% iFs = [5, 6, 7];
load('active_trials.mat', 'Fm')
username = 'timvd';

for iF = 5

Ks = find(Fm(:,iF) > 0.1); % only consider active trials
n = 1; % ISI number
m = 1; % AMP number
tiso = 3; % isometric time (s)

load('parms_v2.mat','sparms','pparms');


Data = get_data(username, iF,n,m,Ks,tiso);
[tis, Cas, Lis, vis, ts] = create_input(tiso, Data.dTt, Data.dTc, Data.ISI, Data.Ca(Ks));
gamma = 108.3333; % length scaling
Liss = Lis * gamma;

figure(6 + iF * 10)
subplot(411)
plot(Data.t, Data.C,'r.'); hold on
plot(tis, Cas, 'b', 'linewidth',1); 
box off

subplot(412)
plot(Data.t, Data.v,'r.'); hold on
plot(tis, vis,'b',  'linewidth',1); 
box off

subplot(413)
plot(Data.t, Data.L,'r.'); hold on
plot(tis, Lis,'b',  'linewidth',1); 
box off

subplot(414)
plot(Data.t, Data.F,'r.'); hold on
box off


for kk = 1:2
    if kk == 1
        parms = sparms(iF);
    else
        parms = pparms(iF);
    end
    
    parms.ti = tis;
    parms.vts = vis;
    parms.Cas = Cas;

    odeopt = odeset('maxstep', 1e-3);
    x0 = 1e-3 * ones(6,1);
    xp0 = zeros(size(x0));

    xsol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 max(parms.ti)], x0, xp0, odeopt);

    xF = (xsol.y(1,:) + xsol.y(2,:)) * parms.Fscale;
    xt = xsol.x;
    xFi = interp1(xt, xF, tis) + parms.Fpe_func(Liss, parms);

    plot(tis, xFi);
end

end