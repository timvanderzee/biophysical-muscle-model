clear all; close all; clc
[username, githubfolder] = get_paths();

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
iFs = [2,3,5,6,7,8,11];
n = [3 1]; % ISI number
m = [7 1]; % AMP number
tiso = 3; % isometric time (s)
N = 500;

versions = {'v3', 'v3'};

mcodes = [1 1 1; 1 2 1];

% specify data
load('active_trials.mat', 'Fm')

SRSrel = nan(2, 100, max(iFs));
SRSreld = nan(100, max(iFs));

for iF = iFs
    
    cd(['C:\Users\',username,'\OneDrive - KU Leuven\9. Short-range stiffness\matlab\data'])
    load([fibers{iF},'_cor_new.mat'],'data')
    
    Ks = find(Fm(:,iF) > .05); % only consider active trials
    Data = prep_data_v2(data,n, m,Ks,tiso);
    [tis, Cas, Lis, vis, ts] = create_input(tiso, Data.dTt, Data.dTc, Data.ISI, Data.Ca(Ks), N);
    
    foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iF}];
    if ~isfolder(foldername)
        mkdir(foldername)
    end
    
    cd(foldername)
    
    figure(1)
    nexttile

    figure(2)
    nexttile
    
    for j = 1:length(versions)
        
        mcode = mcodes(j,:);
        [output_mainfolder, filename, opt_type, ~] = get_folder_and_model(mcode);


        load(['parms_', filename, '_', versions{j}, '.mat'], 'newparms', 'optparms', 'out', 'bnds');
        
        odeopt = odeset('maxstep', 3e-3);
        x0 = 1e-3 * ones(7,1);
        
        if newparms.J1 == 0
            x0(6) = 1;
        end
        
        % simulate
        sol = ode15s(@(t,y) fiber_dynamics_explicit_no_tendon(t,y, newparms), [0 max(newparms.ti)], x0, odeopt);
        
        % extract
        F = (sol.y(1,:) + sol.y(2,:)) * newparms.Fscale;
        t = sol.x;
        Fi = interp1(t, F,  newparms.ti) + newparms.Fpe_func(newparms.Lts, newparms);
        
        % plot
        figure(1)
        plot(Data.t, Data.F, 'k.'); hold on
        plot(newparms.ti, Fi,'-', 'linewidth', 2);
        
        Fii = interp1(newparms.ti, Fi, Data.t);
        Lii = interp1(newparms.ti, newparms.Lts/newparms.gamma, Data.t);
        
        RMSD(iF,j) = mean(abs(Fii - Data.F), 'omitnan');
        Fcost(iF,j) = sum(out.Fcost(:));
        
        % pre-allocate
        ns = nan(length(Ks), 2);
        ds = nan(length(Ks), 2);
        F0 = nan(length(Ks), 1);
        
        [id0,id1,id2] = get_indices(Data.t, tiso, ts, Data.dTt, Data.dTc(1), Data.ISI(1), Data.Ca(Ks));
        
        for i = 1:length(Ks)
            
            np1 = polyfit(Lii(id1(i,:)), Fii(id1(i,:)), 1);
            np2 = polyfit(Lii(id2(i,:)), Fii(id2(i,:)), 1);
            
            dp1 = polyfit(Data.L(id1(i,:)), Data.F(id1(i,:)), 1);
            dp2 = polyfit(Data.L(id2(i,:)), Data.F(id2(i,:)), 1);
            
            ns(i,:) = [np1(1) np2(1)];
            ds(i,:) = [dp1(1) dp2(1)];
            
            F0(i) = mean(Data.F(id0(i,:)));            
 
        end

        %% SRS plot
        figure(2)
        if j == 1
            plot(F0, ds(:,1)./ds(:,2),'ko-'); hold on
        end
        plot(F0, ns(:,1)./ns(:,2),'o-')
        box off
        xlabel('Isometric force (F_0)')
        ylabel('Relative stiffness')
        xlim([0 1.05])
        
        Fs = linspace(0,1,100);
        
        SRSrel(j,:,iF) = interp1(F0, ns(:,1)./ns(:,2), Fs);
        
    end
    
    SRSreld(:,iF) = interp1(F0,ds(:,1)./ds(:,2), Fs);
    
end

figure(2)
legend('Data', 'Model 1', 'Model 2', 'location', 'best')
legend boxoff

%%
figure(10)
plot(Fs, mean(SRSrel(1,:,:), 3, 'omitnan')); hold on
plot(Fs, mean(SRSrel(2,:,:), 3, 'omitnan')); hold on
plot(Fs, mean(SRSreld(:,:), 2, 'omitnan'), 'k-'); hold on