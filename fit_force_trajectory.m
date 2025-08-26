clear all; close all; clc

addpath(genpath('C:\Users\timvd\Documents\muscle-thixotropy'))
addpath(genpath('C:\Users\timvd\Documents\casadi-windows-matlabR2016a-v3.5.5'))

% Import casadi libraries
import casadi.*; 

%% specify data
i = 6;
Ks = 1:4;
n = 3;
m = 7;
tiso = 3;

%% load data
[Data, tis, Lis, vis, Cas, ts, id0, id1, id2] = get_data(i,n,m,Ks,tiso);

% interpolate force
Fis = interp1(Data.t, Data.F, tis);

close all
figure(1)
subplot(411)
plot(Data.t, Data.Ca,'r.'); hold on
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


%% get parameters
mcode = [1 1 1];
vs = {'\', '\'};
% ts = '2_trials';

cd('C:\Users\timvd\Documents\muscle-thixotropy\new_model\get_variable')
[output_mainfolder, filename, opt_type, ~] = get_folder_and_model(mcode);

disp(filename)
output_folder = [opt_type,'\normalized\with_PE_optimized\2_trials'];

output_dir = [output_mainfolder{1}, '\', filename,vs{1}, output_folder];
cd(output_dir)

% get parallel parameters from that fiber
load([filename,'_F', num2str(i),'_best.mat'],'parms','exitflag','fopt','C0','Cbounds','model','P0','P')
C = p_to_c(P, Cbounds);
parms = C_to_parms(C, parms, parms.optvars);
kpe = parms.kpe;
Fpe0 = parms.Fpe0;

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
parms.ti = tis;
parms.vts = vis;
parms.Cas = Cas;
parms.Lts = Lis;

odeopt = odeset('maxstep', 1e-3);
x0 = 1e-3 * ones(6,1);
xp0 = zeros(size(x0));

osol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 max(parms.ti)], x0, xp0, odeopt);
[~,xdot] = deval(osol, osol.x);

gamma = 108.3333; % length scaling
Liss = Lis * gamma;
oF = (osol.y(1,:) + osol.y(2,:)) * parms.Fscale;
ot = osol.x;

oFi = interp1(ot, oF, tis) + parms.Fpe_func(Liss, parms);

subplot(414); hold on
plot(tis, oFi,'b'); hold on

%% intervals of interest
id = nan(length(Ks), length(parms.ti));
% id2 = nan(length(Ks), length(parms.ti));

for k = 1:length(Ks)
    id(k,:) = parms.ti > (ts(k) + 2.2-0.0843) & (parms.ti < ts(k)+3-0.0843);
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
w3 = .1; 	% weight for regularization
w = [w1 w2 w3];

% specify biophysical parameters to be fitted
optparms = {'f', 'k11', 'k22', 'k21', 'JF', 'kon', 'koop', 'kse', 'kse0'};
% optparms = {'f', 'k11', 'k22', 'k21', 'JF', 'kon', 'koop'};

fparms = parms;
fparms.J1 = 50;
fparms.J2 = 200;

if ishandle(3), close(3); end; figure(3)
[newparms, out] = fit_model_parameters_v2(opti, optparms, w, Xdata, fparms);
set(gcf,'units','normalized','position',[.2 .2 .4 .6])

sparms(i) = newparms;

%% visualize fitted parameters
id = i;
clear Y
for j = 1:length(id)
    for i = 1:length(optparms)
        Y(i,j) = eval(['sparms(', num2str(id(j)),').',optparms{i}]);
    end
end

figure(6)
bar(categorical(optparms), Y)
set(gca,'YScale','log')
yline(.1,'r--')
yline(2e3,'r--')

%% test with fitted paramers
parms.ti = tis;
parms.vts = vis;
parms.Cas = Cas;

odeopt = odeset('maxstep', 3e-3);
x0 = 1e-3 * ones(6,1);
xp0 = zeros(size(x0));

nsol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, newparms), [0 max(parms.ti)], x0, xp0, odeopt);
[~,xdot] = deval(nsol, nsol.x);

nF = (nsol.y(1,:) + nsol.y(2,:)) * parms.Fscale;
nt = nsol.x;

nFi = interp1(nt, nF, tis) + parms.Fpe_func(Liss, newparms);

figure(1)
subplot(414); hold on
plot(tis, nFi,'m'); hold on


%% estimate SRS
% lts   = interp1(Data.t, Data.Lf, parms.ti);
nFi = interp1(nsol.x, nF, Data.t)  + parms.Fpe_func(Data.L * gamma, newparms);
oFi = interp1(osol.x, oF, Data.t)  + parms.Fpe_func(Data.L * gamma, parms);

if ishandle(4), close(4); end
figure(4)
subplot(131);
plot(Data.L, Data.F,'.'); hold on

subplot(132);
plot(Data.L, oFi,'.'); hold on

subplot(133)
plot(Data.L, nFi,'.'); hold on

for i = 1:length(Ks)
    
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

%% summary plot
if ishandle(5), close(5); end
figure(5)
plot(F0, ds(:,1)./ds(:,2),'o-'); hold on
plot(F0, os(:,1)./os(:,2),'o-')
plot(F0, ns(:,1)./ns(:,2),'o-')

return
%% test on the condition without conditioning stretch
n = 1;
m = 1;
[Data, tis, Lis, vis, Cas, ts, id0, id1, id2] = get_data(i,n,m,Ks,tiso);

figure(10)
subplot(411)
plot(Data.t, Data.Ca,'r.'); hold on
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

parms.ti = tis;
parms.vts = vis;
parms.Cas = Cas;

dsol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 max(parms.ti)], x0, xp0, odeopt);

newparms.ti = tis;
newparms.vts = vis;
newparms.Cas = Cas;

xsol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, newparms), [0 max(newparms.ti)], x0, xp0, odeopt);

dF = (dsol.y(1,:) + dsol.y(2,:)) * parms.Fscale;
xF = (xsol.y(1,:) + xsol.y(2,:)) * parms.Fscale;

subplot(414); hold on
plot(dsol.x, dF,'b'); hold on
plot(xsol.x, xF,'m'); hold on


%%
function yp=grad5(y,dx)
% function yp=grad5(y,dx)
% 031091 KvS
% purpose: calculation of 5 point derivative
% inputs : y : vector containing signal as a function of x
%          dx: change in x between steps
% output : yp: derivative of y with respect to x
 
if nargin~=2
  disp('OOPS: grad5 expects 2 input arguments')
  disp('      hit any key to continue');pause
  return
end

%y=y(:);
[nrow,ncol]=size(y);
% keyboard

yp(1,1:ncol)=zeros(1,ncol);
yp(2,1:ncol )=zeros(1,ncol);
yp(nrow-1,1:ncol)=zeros(1,ncol);
yp(nrow,1:ncol)=zeros(1,ncol);

yp(1,1:ncol)=(y(2,:)-y(1,:))/dx;
yp(2,1:ncol)=(y(3,:)-y(1,:))/(2*dx);
yp(nrow-1,1:ncol)=(y(nrow,:)-y(nrow-2,:))/(2*dx);
yp(nrow,1:ncol)=(y(nrow,:)-y(nrow-1,:))/dx;

coef=[1 -8 0 8 -1];
for i=3:nrow-2;
  yp(i,:)=(coef*y(i-2:i+2,:))/(12*dx);
end;
end
%if flip; y=y';yp=yp';end;

function[Data, tis, Lis, vis, Cas, ts, id0, id1, id2] = get_data(i,n,m,Ks,tiso)

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

%% process data
% load short-range stiffness data (skinned rat soleus muscle fibers), 
cd('C:\Users\timvd\OneDrive - KU Leuven\9. Short-range stiffness\matlab\data')
load([fibers{i},'_cor_new.mat'],'data')

Data.F = [];
Data.L = [];
Data.t = [];
Data.Ca = [];

dTt = .0383/.4545; % test stretch (= constant)

for k = 1:length(Ks)

    F = data.Fexp(:,Ks(k),n,m);
    L = data.Lexp(:,Ks(k),n,m);
    t = data.texp(:,Ks(k),n,m);
    Ca = 10^(6-data.pCas(Ks(k))) * ones(size(F));

    ISI = data.ISIs(n)/1000;
    dTc = data.AMPs(m)/10000 / .4545; % conditioning stretch
    
    trange = [-.8 .2];
    id1 = t > trange(1) & t < trange(2);
   
    % save
    Data.t = [Data.t; t(id1) + k * tiso - 3*dTt - .005]; % move to agree with idealized input
    Data.F = [Data.F; F(id1)];
    Data.L = [Data.L; L(id1)];
    Data.Ca = [Data.Ca; Ca(id1)];
       
end

% intervals
ts = 0:tiso:(tiso*(length(Ks)-1));

%% get velocity (through interpolating and filtering)
ti = linspace(0,max(Data.t), 10000);
Li = interp1([0; Data.t], [0; Data.L], ti);

% filter
dt = mean(diff(ti),'omitnan');
fs = 1/dt;
Wn = 30 / (.5 * fs);
[b,a] = butter(2, Wn);

Lf = Li;
idf = isfinite(Li);

Lf(idf) = filtfilt(b,a, Li(idf));

vf = grad5(Lf(:), dt);

Data.v = interp1(ti, vf, Data.t);

%% create idealized input
Ts = cumsum([tiso-(dTt*3+dTc*2+ISI); dTc; dTc; ISI; dTt; dTt; dTt]);
tx = linspace(0, Ts(end), 5000);

vx = nan(size(tx));
vx(tx <= Ts(1)) = 0;
vx(tx > Ts(1) & tx <=  Ts(2)) = .4545;
vx(tx > Ts(2) & tx <=  Ts(3)) = -.4545;
vx(tx > Ts(3) & tx <=  Ts(4)) = 0;
vx(tx > Ts(4) & tx <=  Ts(5)) = .4545;
vx(tx > Ts(5) & tx <=  Ts(6)) = 0;
vx(tx > Ts(6) & tx <=  Ts(7)) = -.4545;

Lx = cumtrapz(tx, vx);

ti = linspace(0, Ts(end), 500);
vi = interp1(tx, vx, ti);
Li = interp1(tx, Lx, ti);
%%
Ca = 10.^(6-data.pCas);

vis = [];
Cas = [];
Lis = [];

for k = 1:length(Ks)
    vis = [vis vi];
    Cas = [Cas Ca(Ks(k)) * ones(size(vi))];
    Lis = [Lis Li];
end

tis = linspace(0, Ts(end)*length(Ks), length(ti)*length(Ks));
% Lis = cumtrapz(tis, vis);

%% SRS indices
id0 = nan(length(Ks), 100);
id1 = nan(length(Ks), 10);
id2 = nan(length(Ks), 10);

for i = 1:length(Ks)
    id0(i,:) = find(Data.t < (ts(i) + 3 - 3*dTt - 2*dTc - ISI),100, 'last');
    id1(i,:) = find(Data.t > (ts(i) + 3 - 3*dTt),10, 'first');
    id2(i,:) = find(Data.t > (ts(i) + 3 - 3*dTt - 2*dTc - ISI),10, 'first');
end

%%


end
