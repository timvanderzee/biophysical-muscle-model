clear all; close all; clc
% mainfolder = 'C:\Users\timvd\Documents';
mainfolder = 'C:\Users\u0167448\Documents\';
addpath(genpath([mainfolder, 'GitHub\muscle-thixotropy']))
addpath(genpath([mainfolder, 'GitHub\biophysical-muscle-model']))

% fibers
iFs = [1 2 3, 5, 6, 7, 8, 10, 11];
% iFs = 6;

% get parameters
load('parms_v4.mat','sparms','pparms');

x0 = 1e-3 * ones(6,1);
xp0 = zeros(size(x0));
odeopt = odeset('maxstep', 5e-2);
gamma = 108.3333; % length scaling

tiso = 5;

%% step 1: force - pCa
pCas = flip([9, 7:-.2:5, 4.5]);

ks = 1;

Fss = nan(length(pCas), 2, iFs(end));

for kk = ks
    if kk == 1
        aparms = sparms;
    else
        aparms = pparms;
    end
     
    for iF = iFs
        
        parms = aparms(iF);
        
        Ca = 10.^(-pCas+6);
        
        [tis, Cas, Lis, vis, ts] = create_input(tiso, 0, 0, 0, Ca, 2000);
        id0 = get_indices(tis, tiso, ts, 0, 0, 0, Ca);
        
        parms.ti = tis;
        parms.vts = vis;
        parms.Cas = Cas;
        parms.Lts = Lis;
        
        tic
        sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 max(tis)], x0, xp0, odeopt);
        toc
        
        oF = (sol.y(1,:) + sol.y(2,:)) * parms.Fscale;
        ot = sol.x;
        
        oFi = interp1(ot, oF, tis);
        
        figure(1)
        subplot(411)
        plot(tis, Cas, 'b', 'linewidth',1);
        
        subplot(412)
        plot(tis, vis,'b',  'linewidth',1);
        
        subplot(413)
        plot(tis, Lis,'b',  'linewidth',1);
        
        subplot(414);
        plot(tis, oFi, tis(id0), oFi(id0),'r.'); hold on
        
        for j = 1:size(id0,1)
            Fss(j,kk,iF) = mean(oFi(id0(j,:)));
        end
    end
end

%% find pCa to yield specified force
if ishandle(2), close(2); end

Fsi = linspace(.99,.01,7);
pCai = nan(2, length(Fsi), iFs(end));

for iF = iFs
    figure(2)
    nexttile
    for kk = ks
        pCai(kk,:,iF) = interp1(Fss(:,kk,iF), pCas, Fsi, [], 'extrap');
        
        plot(pCas, Fss(:,kk,iF), pCai(kk,:,iF), Fsi, 'o'); hold on
    end
end


%% step 2: stretch-shortening
% iFs = [5, 6, 7];

% conditions
ISIs = logspace(-2,0,5);
AMPs = linspace(0.01, .0383, 5);

% outputs
F0      = nan(length(Fsi), length(ISIs),length(AMPs), 11, 2);
Scond   = nan(length(Fsi), length(ISIs),length(AMPs), 11, 2);
Stest   = nan(length(Fsi), length(ISIs),length(AMPs), 11, 2);

for kk = ks
    
    if kk == 1
        aparms = sparms;
    else
        aparms = pparms;
    end
    
    for iF = iFs
        pCas = pCai(kk,:,iF);
        
        % avoid issues
        pCas(pCas < 4) = 4;
        
        parms = aparms(iF);
        
        for ii = 1:length(AMPs)
            
            AMP = AMPs(ii);
            dTt = .0383/.4545; % test stretch (= constant)
            dTc = AMP / .4545; % conditioning stretch
            
            for jj = 1:length(ISIs)
                
         
                ISI = ISIs(jj);
                Ca = 10.^(-pCas+6);
                
                [tis, Cas, Lis, vis, ts] = create_input(tiso, dTt, dTc, ISI, Ca, 2000);
                
                [id0,id1,id2] = get_indices(tis, tiso, ts, dTt, dTc, ISI, Ca);
                
                parms.ti = tis;
                parms.vts = vis;
                parms.Cas = Cas;
                parms.Lts = Lis;

                figure(3)
                subplot(411)
                plot(tis, Cas, 'b', 'linewidth',1);
                
                subplot(412)
                plot(tis, vis,'b',  'linewidth',1);
                
                subplot(413)
                plot(tis, Lis,'b',  'linewidth',1);
                
                % run simulation
                tic
                sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 max(tis)], x0, xp0, odeopt);
                toc
                
                Liss = Lis * gamma;
                oF = (sol.y(1,:) + sol.y(2,:)) * parms.Fscale;
                ot = sol.x;
                
                oFi = interp1(ot, oF, tis) + parms.Fpe_func(Liss, parms);
                
                subplot(414);
                plot(tis, oFi,'b', tis(id0), oFi(id0),'r.',tis(id1), oFi(id1),'k.', tis(id2), oFi(id2),'g.')
                
                for j = 1:4
                    subplot(4,1,j)
                    box off
                    
                end
                
                %
                for i = 1:length(Ca)
                    
                    np1 = polyfit(Lis(id1(i,:)), oFi(id1(i,:)), 1);
                    np2 = polyfit(Lis(id2(i,:)), oFi(id2(i,:)), 1);
                    
                    Stest(i,jj,ii,iF,kk) = np1(1);
                    Scond(i,jj,ii,iF,kk) = np2(1);
                    
                    F0(i,jj,ii,iF,kk) = mean(oFi(id0(i,:)));
                    
                end
            end
        end
    end
end

