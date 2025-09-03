clear all; close all; clc

mainfolder = 'C:\Users\u0167448\Documents\GitHub';
% mainfolder = 'C:\Users\timvd\Documents';
% username = 'timvd';
username = 'u0167448';

addpath(genpath([mainfolder, '\muscle-thixotropy']))
addpath(genpath([mainfolder, '\casadi-windows-matlabR2016a-v3.5.5']))
addpath(genpath([mainfolder, '\biophysical-muscle-model']))

% Import casadi libraries
import casadi.*; 

%% specify data
load('active_trials.mat', 'Fm')
iFs = [1 2 3, 5, 6, 7, 8, 10, 11];
% iFs = [5, 6, 7];

for iF = iFs

% iF = 6; % fiber number (note: #4 and #9 lack pCa = 4.5)
Ks = find(Fm(:,iF) > 0.1); % only consider active trials
n = 3; % ISI number
m = 7; % AMP number
tiso = 3; % isometric time (s)

%% load data
Data = get_data(username, iF,n,m,Ks,tiso);
[tis, Cas, Lis, vis, ts] = create_input(tiso, Data.dTt, Data.dTc, Data.ISI, Data.Ca(Ks));

% interpolate force
Fis = interp1(Data.t, Data.F, tis);

figure(1 + iF*10)
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
plot(Data.t, Data.F,'r.');
box off

for j = 1:4
     subplot(4,1,j)
     box off
    
    for i = 1:length(ts)
        xline(ts(i),'k--')
        xline(ts(i) + mean(diff(ts)) - 3*Data.dTt - 2*Data.dTc - Data.ISI, 'k:')
        xline(ts(i) + mean(diff(ts)) - 3*Data.dTt, 'k:')
        
    end
end

%% get parameters
mcode = [1 1 1];
vs = {'\', '\'};

cd([mainfolder, '\muscle-thixotropy\new_model\get_variable'])
[output_mainfolder, filename, opt_type, ~] = get_folder_and_model(mcode);

disp(filename)
output_folder = [opt_type,'\normalized\with_PE_optimized\2_trials'];

output_dir = [output_mainfolder{1}, '\', filename,vs{1}, output_folder];
cd(output_dir)

% get parallel parameters from that fiber
load([filename,'_F', num2str(iF),'_best.mat'],'parms','exitflag','fopt','C0','Cbounds','model','P0','P')
C = p_to_c(P, Cbounds);
oparms = parms;
oparms = C_to_parms(C, oparms, oparms.optvars);
kpe = oparms.kpe;
Fpe0 = oparms.Fpe0;
oparms.act = 1;
oparms.Noverlap = 1;
oparms.Fscale = 1.5;

% get other parameters from other fiber
ii = 6; % fiber from which parameters are obtained
load([filename,'_F', num2str(ii),'_best.mat'],'parms','exitflag','fopt','C0','Cbounds','model','P0','P')
C = p_to_c(P, Cbounds);
parms = C_to_parms(C, parms, parms.optvars);
parms = calc_dependent_parms(parms);  

parms.act = 1;
parms.Noverlap = 1;
parms.Fscale = 1.5;
parms.kpe = kpe;
parms.Fpe0 = Fpe0;

%% evaluate
% note: not required for fitting
oparms.ti = tis;
oparms.vts = vis;
oparms.Cas = Cas;
oparms.Lts = Lis;

odeopt = odeset('maxstep', 1e-3);
x0 = 1e-3 * ones(6,1);
xp0 = zeros(size(x0));

osol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, oparms), [0 max(oparms.ti)], x0, xp0, odeopt);

gamma = 108.3333; % length scaling
Liss = Lis * gamma;
oF = (osol.y(1,:) + osol.y(2,:)) * oparms.Fscale;
ot = osol.x;

oFi = interp1(ot, oF, tis) + oparms.Fpe_func(Liss, parms);

subplot(414); hold on
plot(tis, oFi,'b'); hold on

%% intervals of interest
id = nan(length(Ks), length(oparms.ti));

for k = 1:length(Ks)
    id(k,:) = oparms.ti > (ts(k) + 2.2-0.0843) & (oparms.ti < ts(k)+3-0.0843);
end

idF = find(sum(id,1) & isfinite(Fis));

subplot(414); hold on
plot(tis(idF), oFi(idF),'g.'); hold on

%% select data for fitting
Xdata.t = tis;
Xdata.F = Fis;
Xdata.v = vis;
Xdata.L = Liss;
Xdata.Cas = Cas;
Xdata.idF = idF;
Xdata.idC = idF;

%% do fitting
% initialise opti structure
opti = casadi.Opti(); 

% define weigth vector
w1 = 100;    % weight for fitting force-velocity
w2 = 100;   % weight for fitting short-range stiffness
w3 = 1; 	% weight for regularization
w = [w1 w2 w3];

% specify biophysical parameters to be fitted
optparms = {'f', 'k11', 'k22', 'k21', 'JF', 'J1','J2', 'kon', 'koop', 'kse', 'kse0'};

fparms = parms;
% fparms.J1 = 50;
% fparms.J2 = 200;

figure(2 + iF*10)
[newparms, out] = fit_model_parameters_v2(opti, optparms, w, Xdata, fparms);
set(gcf,'units','normalized','position',[.2 .2 .4 .6])

sparms(iF) = newparms;
pparms(iF) = oparms;

end

%% visualize fitted parameters
for i = 1:length(optparms)
    Y(i,1) = eval(['pparms(', num2str(iF),').',optparms{i}]);
    Y(i,2) = eval(['fparms.',optparms{i}]);
    Y(i,3) = eval(['sparms(', num2str(iF),').',optparms{i}]);
end

figure(1)
nexttile 
bar(categorical(optparms), Y)
set(gca,'YScale','log')
yline(.1,'r--')
yline(2e3,'r--')
legend('Old','IG','New','location','best')

%% test with fitted paramers
Kss = [Ks; 7]; % only consider active trials
Data = get_data(username, iF,n,m,Kss,tiso);
[tis, Cas, Lis, vis, ts] = create_input(tiso, Data.dTt, Data.dTc, Data.ISI, Data.Ca(Kss));
Liss = Lis * gamma;

newparms.ti = tis;
newparms.vts = vis;
newparms.Cas = Cas;

odeopt = odeset('maxstep', 3e-3);
x0 = 1e-3 * ones(6,1);
xp0 = zeros(size(x0));

nsol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, newparms), [0 max(tis)], x0, xp0, odeopt);
[~,xdot] = deval(nsol, nsol.x);

nF = (nsol.y(1,:) + nsol.y(2,:)) * parms.Fscale;
nt = nsol.x;

nFi = interp1(nt, nF, tis) + parms.Fpe_func(Liss, newparms);

figure(1 + iF * 10)
subplot(414); hold on
plot(tis, nFi,'m'); hold on

%% estimate SRS
% lts   = interp1(Data.t, Data.Lf, parms.ti);
nFi = interp1(nsol.x, nF, Data.t)  + parms.Fpe_func(Data.L * gamma, newparms);
oFi = interp1(osol.x, oF, Data.t)  + parms.Fpe_func(Data.L * gamma, parms);

figure(4 + iF * 10)
subplot(131);
plot(Data.L, Data.F,'.'); hold on

subplot(132);
plot(Data.L, oFi,'.'); hold on

subplot(133)
plot(Data.L, nFi,'.'); hold on

% pre-allocate
ns = nan(length(Kss), 2);
os = nan(length(Kss), 2);
ds = nan(length(Kss), 2);
F0 = nan(length(Kss), 1);

[id0,id1,id2] = get_indices(Data.t, ts, Data.dTt, Data.dTc, Data.ISI, Data.Ca(Kss));

for i = 1:length(Kss)
    
    np1 = polyfit(Data.L(id1(i,:)), nFi(id1(i,:)), 1);
    np2 = polyfit(Data.L(id2(i,:)), nFi(id2(i,:)), 1);
    
    op1 = polyfit(Data.L(id1(i,:)), oFi(id1(i,:)), 1);
    op2 = polyfit(Data.L(id2(i,:)), oFi(id2(i,:)), 1);
    
    dp1 = polyfit(Data.L(id1(i,:)), Data.F(id1(i,:)), 1);
    dp2 = polyfit(Data.L(id2(i,:)), Data.F(id2(i,:)), 1);
    
    ns(i,:) = [np1(1) np2(1)];
    os(i,:) = [op1(1) op2(1)];
    ds(i,:) = [dp1(1) dp2(1)];
    
    F0(i) = mean(Data.F(id0(i,:)));
    
    subplot(131)
    plot(Data.L(id0(i,:)), Data.F(id0(i,:)), 'g.')
    plot(Data.L(id1(i,:)), Data.F(id1(i,:)), 'r.')
    plot(Data.L(id2(i,:)), Data.F(id2(i,:)), 'y.')
    
    subplot(132)
    plot(Data.L(id1(i,:)), oFi(id1(i,:)), 'r.')
    plot(Data.L(id2(i,:)), oFi(id2(i,:)), 'y.')
    
    subplot(133)
    plot(Data.L(id1(i,:)), nFi(id1(i,:)), 'r.')
    plot(Data.L(id2(i,:)), nFi(id2(i,:)), 'y.')
end

titles = {'Data','Old parameters', 'New parameters'};
for i = 1:3
    subplot(1,3,i)
    box off
    xlabel('Length (L_0)')
    ylabel('Force (F_0)')
    title(titles{i})
end

%% summary plot
figure(2)
nexttile
plot(F0, ds(:,1)./ds(:,2),'o-'); hold on
plot(F0, os(:,1)./os(:,2),'o-')
plot(F0, ns(:,1)./ns(:,2),'o-')
box off
xlabel('Isometric force (F_0)')
ylabel('Relative stiffness')
legend('Data', 'Old parameters', 'New parameters', 'location', 'best')
legend boxoff
xlim([0 1.05])

%% test on the condition without conditioning stretch
% n = 1;
% m = 1;
% [Data, tis, Lis, vis, Cas, ts, id0, id1, id2] = get_data(iF,n,m,Ks,tiso);
% 
% figure(6 + iF * 10)
% subplot(411)
% plot(Data.t, Data.Ca,'r.'); hold on
% plot(tis, Cas, 'b', 'linewidth',1); 
% box off
% 
% subplot(412)
% plot(Data.t, Data.v,'r.'); hold on
% plot(tis, vis,'b',  'linewidth',1); 
% box off
% 
% subplot(413)
% plot(Data.t, Data.L,'r.'); hold on
% plot(tis, Lis,'b',  'linewidth',1); 
% box off
% 
% subplot(414)
% plot(Data.t, Data.F,'r.'); hold on
% box off
% 
% parms.ti = tis;
% parms.vts = vis;
% parms.Cas = Cas;
% 
% osol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 max(parms.ti)], x0, xp0, odeopt);
% 
% newparms.ti = tis;
% newparms.vts = vis;
% newparms.Cas = Cas;
% 
% oF = (osol.y(1,:) + osol.y(2,:)) * parms.Fscale;
% ot = osol.x;
% oFi = interp1(ot, oF, tis) + parms.Fpe_func(Liss, parms);
% 
% xsol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, newparms), [0 max(newparms.ti)], x0, xp0, odeopt);
% 
% xF = (xsol.y(1,:) + xsol.y(2,:)) * parms.Fscale;
% xt = xsol.x;
% xFi = interp1(xt, xF, tis) + parms.Fpe_func(Liss, newparms);
% 
% 
% plot(tis, oFi,'b'); hold on
% plot(tis, xFi,'m');
% end

% %% compare parameter values
% figure(3)
% 
% Ys = nan(length(optparms), iFs(end));
% Yp = nan(length(optparms), iFs(end));
% 
% for iF = iFs
%     for i = 1:length(optparms)
%         Ys(i,iF) = eval(['sparms(', num2str(iF),').',optparms{i}]);
%         Yp(i,iF) = eval(['pparms(', num2str(iF),').',optparms{i}]);
%     end
% end

%%

%if flip; y=y';yp=yp';end;

