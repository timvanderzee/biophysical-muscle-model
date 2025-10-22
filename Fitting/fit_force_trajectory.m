clear all; close all; clc

usernames = {'timvd','u0167448'};

for i = 1:length(usernames)
    docfolder =  ['C:\Users\', usernames{i}, '\Documents'];
    
    if isfolder(docfolder)
        mainfolder = docfolder;
        username = usernames{i};
    end
end

if contains(mainfolder, 'timvd')
    githubfolder = mainfolder;
else
    githubfolder = [mainfolder, '\GitHub'];
end

casadifolder = 'C:\Users\u0167448\Documents\casadi-3.7.1-windows64-matlab2018b';

addpath(genpath([githubfolder, '\muscle-thixotropy\new_model\']))
addpath(genpath(casadifolder))
% rmpath(genpath('C:\Users\u0167448\Documents\casadi-3.7.1-windows64-matlab2018b'))
% addpath(genpath('C:\GBW_MyPrograms\casadi-3.6.7-windows64-matlab2018b'))
addpath(genpath([githubfolder, '\biophysical-muscle-model']))

%% specify data
load('active_trials.mat', 'Fm')
% iFs = [1 2 3, 5, 6, 7, 8, 10, 11];
iFs = [2,3,5,6,7,8,11];
iFs = 6;

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

for iF = iFs

% iF = 6; % fiber number (note: #4 and #9 lack pCa = 4.5)
Ks = find(Fm(:,iF) > .05); % only consider active trials
% Ks = find(isfinite(Fm(:,iF))); % only consider active trials
n = [3 1]; % ISI number
m = [7 1]; % AMP number
tiso = 3; % isometric time (s)

%% load data
cd(['C:\Users\',username,'\OneDrive - KU Leuven\9. Short-range stiffness\matlab\data'])
load([fibers{iF},'_cor_new.mat'],'data')

Data = prep_data_v2(data,n, m,Ks,tiso);
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
plot(Data.t, Data.F,'r.'); hold on
% plot(tis, Fis, 'r--')
box off

for j = 1:4
     subplot(4,1,j)
     box off
    
%     for i = 1:length(ts)
%         xline(ts(i),'k--')
%         xline(ts(i) + mean(diff(ts)) - 3*Data.dTt - 2*Data.dTc - Data.ISI, 'k:')
%         xline(ts(i) + mean(diff(ts)) - 3*Data.dTt, 'k:')
%     end
end

%% get parameters
mcode = [1 1 1];
vs = {'\', '\'};
% vs = {'\full\'};

cd([githubfolder, '\muscle-thixotropy\new_model\get_variable'])
[output_mainfolder, filename, opt_type, ~] = get_folder_and_model(mcode);

disp(filename)
output_folder = [opt_type,'\normalized\with_PE_optimized\2_trials'];

output_dir = [output_mainfolder{1}, '\', filename,vs{1}, output_folder];
cd(output_dir)

% get other parameters from other fiber
% ii = 6; % fiber from which parameters are obtained
load([filename,'_F', num2str(iF),'_best.mat'],'parms','exitflag','fopt','C0','Cbounds','model','P0','P')
C = p_to_c(P, Cbounds);
parms = C_to_parms(C, parms, parms.optvars);
parms = calc_dependent_parms(parms);  

parms.act = 1;
parms.Noverlap = 1;
parms.approx = 0;
parms.K = 100;
oparms = parms;

%% evaluate current values
h = 10e-9; % powerstroke size
s = 2.6e-6; % sarcomere length

% set some parameters
parms.ti = tis;
parms.vts = vis;
parms.Cas = Cas;
parms.Lts = Lis;
parms.gamma = .5*s / h; % length scaling
parms.J1 = 6.17;
parms.koop = 5.7;
parms.JF = 1e3;

% get initial guess
IG = get_initial_guess(tis, Cas, vis, parms);

Liss = Lis * parms.gamma;
% oF = (osol.y(1,:) + osol.y(2,:)) * parms.Fscale;
% ot = osol.x;

oFi = (IG.Q0i+IG.Q1i) * parms.Fscale + parms.Fpe_func(Liss, parms);

subplot(414); hold on
plot(tis, oFi,'b'); hold on

%% intervals of interest
id1 = nan(length(Ks), length(parms.ti), length(m));
id2 = nan(length(Ks), length(parms.ti), length(m));
tmax = [0 tiso * length(Ks)];

for i = 1:length(m)
    for k = 1:length(Ks)
        id1(k,:,i) = (parms.ti > (tmax(i) + ts(k) + tiso - 4*Data.dTt - 2*Data.dTc(i) - Data.ISI(i))) & (parms.ti < (tmax(i) + ts(k)+tiso - Data.dTt));
        
        id2(k,:,i) = (parms.ti > (tmax(i) + ts(k) + tiso - 5*Data.dTt - 2*Data.dTc(i) - Data.ISI(i))) & (parms.ti < (tmax(i) + ts(k) + tiso - 3.5*Data.dTt - 2*Data.dTc(i) - Data.ISI(i)));
    end
end

ID1 = sum(id1,3, 'omitnan');
ID2 = sum(id2,3, 'omitnan');

idF = find(sum(ID1,1) & isfinite(Fis));
idC = find(sum(ID2,1));

subplot(414); hold on
plot(tis(idF), oFi(idF),'bx'); hold on
plot(tis(idC), oFi(idC),'rx'); hold on

%% select data for fitting
Xdata.t = tis;
Xdata.F = Fis;
Xdata.v = vis;
Xdata.L = Liss;
Xdata.Cas = Cas;
Xdata.idF = idF;
Xdata.idC = idC;

%% do fitting
bnds.f = [1 200];
bnds.k11 = [1e-2 200];
bnds.k22 = [0 1];
bnds.k21 = [1 200];
bnds.kF = [1 1e4];
bnds.J1 = [1e-3 200];
bnds.J2 = [1 1e3];
bnds.kon = [1 200];
bnds.kse = [1e-3 1];
bnds.kse0 = [1e-4 1];
bnds.koop = [1 200];

% initialise opti structure
opti = casadi.Opti(); 

% define weigth vector
w1 = 100;    % weight for fitting force-velocity
w2 = 100;   % weight for fitting short-range stiffness
w3 = 1000; 	% weight for regularization
w = [w1 w2 w3];

% specify biophysical parameters to be fitted
optparms = {'f', 'k11', 'k22', 'k21', 'J2', 'kon', 'kse', 'kse0'};

parms.kF = parms.J1 * parms.JF;
fparms = parms;

% fparms.k = 1000;

figure(2 + iF*10)
[newparms, out] = fit_model_parameters_v2(opti, optparms, w, Xdata, fparms, IG, bnds);
set(gcf,'units','normalized','position',[.2 .2 .4 .6])

sparms(iF) = newparms;
pparms(iF) = parms;

%% Test the result: run a forward simulation (sanity check)
% [t,x] = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 max(toc)], sol0.y(:,end), xp0, odeopt);
% F = (x(:,1) + x(:,2)) * parms.Fscale;
% Fn = interp1(t, F, toc) + parms.kpe * Lts + parms.Fpe0;

%% visualize fitted parameters

% for i = 1:length(optparms)
%     Y(i,1) = eval(['pparms(', num2str(iF),').',optparms{i}]);
%     Y(i,2) = eval(['fparms.',optparms{i}]);
%     Y(i,3) = eval(['sparms(', num2str(iF),').',optparms{i}]);
% end

clear Y
for i = 1:length(optparms)
    Y(i,1) = (newparms.(optparms{i}) - bnds.(optparms{i})(1)) / diff(bnds.(optparms{i}));
end

figure(1)
% nexttile 
bar(categorical(optparms), Y)
% set(gca,'YScale','log')
% ylrine(.1,'r--')
% yline(2e3,'r--')
% legend('Old','IG','New','location','best')

%% test with fitted paramers
% n = [3 1]; % ISI number
% m = [7 1]; % AMP number

Kss = [Ks]; % only consider active trials
% Data = prep_data(username, iF,n,m,Kss,tiso);
Data = prep_data_v2(data,n,m,Kss,tiso);

[tis, Cas, Lis, vis, ts] = create_input(tiso, Data.dTt, Data.dTc, Data.ISI, Data.Ca(Kss));
Liss = Lis * parms.gamma;

oparms.ti = tis;
oparms.vts = vis;
oparms.Cas = Cas;

newparms.ti = tis;
newparms.vts = vis;
newparms.Cas = Cas;

odeopt = odeset('maxstep', 3e-3);
x0 = 1e-3 * ones(7,1);
xp0 = zeros(size(x0));

oparms.gamma = 108.3;
% parms.approx = 0;
osol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, oparms), [0 max(tis)], x0, xp0, odeopt);

newparms.approx = 1;
newparms.JF = newparms.kF / newparms.J1;
nsol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, newparms), [0 max(tis)], x0, xp0, odeopt);
% [~,xdot] = deval(nsol, nsol.x);

nF = (nsol.y(1,:) + nsol.y(2,:)) * parms.Fscale;
oF = (osol.y(1,:) + osol.y(2,:)) * oparms.Fscale;

nt = nsol.x;
ot = osol.x;

nFi = interp1(nt, nF, tis) + parms.Fpe_func(Liss, newparms);
oFi = interp1(ot, oF, tis) + parms.Fpe_func(Liss, newparms);

figure(1 + iF * 10)
subplot(414); hold on
plot(tis, nFi,'m'); hold on
plot(tis, oFi,'g'); hold on


%%
% close all

% figure(100)
% color = get(gca,'colororder');
% 
% subplot(211)
% plot(Data.t, Data.F, 'k.'); hold on
% plot(nsol.x, nF, 'color', color(1,:)); hold on
% plot(osol.x, oF, 'color', color(2,:));
% 
% oDF = (Data.F - interp1(osol.x, oF, Data.t)).^2;
% nDF = (Data.F - interp1(nsol.x, nF, Data.t)).^2;
% 
% subplot(212)
% plot(Data.t, [nDF(:) oDF(:)],'.')
% 
% %% compare ripped
% if ishandle(100), close(100); end; figure(100)
% subplot(211)
% plot(osol.x, oF); hold on
% plot(nsol.x, nF)
% 
% subplot(212)
% plot(osol.x, osol.y(end,:)); hold on
% plot(nsol.x, nsol.y(end,:))

%% estimate SRS
% lts   = interp1(Data.t, Data.Lf, parms.ti);
nFi = interp1(nsol.x, nF, Data.t)  + parms.Fpe_func(Data.L * parms.gamma, newparms);
oFi = interp1(osol.x, oF, Data.t)  + parms.Fpe_func(Data.L * parms.gamma, parms);

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

[id0,id1,id2] = get_indices(Data.t, tiso, ts, Data.dTt, Data.dTc(1), Data.ISI(1), Data.Ca(Kss));

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
% Data = prep_data_v2(data,n,m,Ks,tiso);
% [tis, Cas, Lis, vis, ts] = create_input(tiso, Data.dTt, Data.dTc, Data.ISI, Data.Ca(Ks));
% Liss = Lis * parms.gamma;
% 
% figure(6 + iF * 10)
% subplot(411)
% plot(Data.t, Data.C,'r.'); hold on
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

%%
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
% subplot(414)
% plot(tis, oFi,'b'); hold on
% plot(tis, xFi,'m');



%% save
% foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iF}];
% if ~isfolder(foldername)
%     mkdir(foldername)
% end
% 
% cd(foldername)
% save(['parms_', filename, '.mat'], 'newparms')

end

%%

%if flip; y=y';yp=yp';end;

