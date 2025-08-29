clear all; close all; clc
mainfolder = 'C:\Users\timvd\Documents';
addpath(genpath([mainfolder, '\muscle-thixotropy']))
addpath(genpath([mainfolder, '\casadi-windows-matlabR2016a-v3.5.5']))
addpath(genpath([mainfolder, '\biophysical-muscle-model']))


ISIs = logspace(-2,0,10);
pCas = flip([9, 7:-.1:6, 4.5]);
pCas = [6.4 6];


AMPs = [.02 .0383];

F0 = nan(length(pCas), length(ISIs),length(AMPs), 11, 2);
Scond = nan(length(pCas), length(ISIs),length(AMPs), 11, 2);
Stest = nan(length(pCas), length(ISIs),length(AMPs), 11, 2);

iFs = [1 2 3, 5, 6, 7, 8, 10, 11];
iFs = 6;


for ii = 1:length(AMPs)

AMP = AMPs(ii);
dTt = .0383/.4545; % test stretch (= constant)
dTc = AMP / .4545; % conditioning stretch



for jj = 1:length(ISIs)
    
tiso = 3;
ISI = ISIs(jj);
Ca = 10.^(-pCas+6);

[tis, Cas, Lis, vis, ts] = create_input(tiso, dTt, dTc, ISI, Ca, 2000);

[id0,id1,id2] = get_indices(tis, tiso, ts, dTt, dTc, ISI, Ca);


%% get parameters
load('parms_v2.mat','sparms','pparms');


for iF = iFs
    
for kk = 1
    if kk == 1
        parms = sparms(iF);
    else
        parms = pparms(iF);
    end

figure(kk)
subplot(411)
plot(tis, Cas, 'b', 'linewidth',1); 

subplot(412)
plot(tis, vis,'b',  'linewidth',1); 

subplot(413)
plot(tis, Lis,'b',  'linewidth',1); hold on

%% run simulation
x0 = 1e-3 * ones(6,1);
xp0 = zeros(size(x0));

% isometric
tic
parms.vts = [0 0];
parms.ti = [0 100];
parms.Cas = Cas(1) * [1 1];
sol0 = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 max(parms.ti)], x0, xp0, []);

parms.ti = tis;
parms.vts = vis;
parms.Cas = Cas;
parms.Lts = Lis;
odeopt = odeset('maxstep', 5e-2);

sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 max(tis)], sol0.y(:,end), xp0, odeopt);
toc

gamma = 108.3333; % length scaling
Liss = Lis * gamma;
oF = (sol.y(1,:) + sol.y(2,:)) * parms.Fscale;
ot = sol.x;

oFi = interp1(ot, oF, tis); %+ parms.Fpe_func(Liss, parms);

subplot(414); hold on
plot(tis, oFi,'b'); hold on
plot(tis(id0), oFi(id0),'r.',  'linewidth',1); 
plot(tis(id1), oFi(id1),'k.',  'linewidth',1); 
plot(tis(id2), oFi(id2),'g.',  'linewidth',1); 

for j = 1:4
     subplot(4,1,j)
     box off
    
    for i = 1:length(ts)
        xline(ts(i),'k--')
    end
end

%%
for i = 1:length(Ca)
    
    np1 = polyfit(Lis(id1(i,:)), oFi(id1(i,:)), 1);
    np2 = polyfit(Lis(id2(i,:)), oFi(id2(i,:)), 1);
    
    Stest(i,jj,ii,iF,kk) = np1(1);
    Scond(i,jj,ii,iF,kk) = np2(1);
    
%     ns = [np1(1) np2(1)];

    F0(i,jj,ii,iF,kk) = mean(oFi(id0(i,:)));

end
end
end
end
end

%% effect of ISI
SRSrel = Stest./Scond(:,:,2,:,:);
close all
figure(10)
color = get(gca,'colororder');

for iF = iFs

    subplot(131);
    plot(squeeze(F0(:,1,:,iF, kk)), squeeze(SRSrel(:,1,:,iF,kk)), 'color', color(kk,:)); hold on
    
    subplot(132)
    semilogx(ISIs, squeeze(SRSrel(1,:,:,iF,kk)), 'color', color(kk,:)); hold on
   
    subplot(133)  
    plot(AMPs, squeeze(SRSrel(:,1,:,iF,kk))', 'color', color(kk,:)); hold on
   
end


return

%%
%% summary figure

% figure(10)
% color = get(gca,'colororder');
% 
% for iF = iFs
%     for kk = 1
% 
%         figure(10)
%         plot(F0(:,1,2,iF, kk), SRSrel(:,1,2,iF,kk), 'color', color(kk,:)); hold on
%     end
% end


%% interpolate
% F0s = linspace(0,1,13);
% SRSi = nan(size(SRSrel));
% 
% for iF = iFs
%     for kk = 1
%         SRSi(:,1,2,iF,kk) = interp1(F0(:,1,2,iF,kk), SRSrel(:,1,2,iF,kk), F0s(:), [], 'extrap');
%     end
% end
% 
% figure(11)
% for kk = 1:2
% plot(F0s, mean(SRSi(:,1,2,:,kk), 2,'omitnan'),'color',color(kk,:)); hold on
% end
% 





