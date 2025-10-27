clear all; close all; clc
[username, githubfolder] = get_paths();
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

% model to be fitted
mcode = [1 2 1];

% settings
save_results = 1;
iFs = 11; %[2,3,5,6,7,8,11];
n = [3 1]; % ISI number
m = [7 1]; % AMP number
tiso = 3; % isometric time (s)

% bounds
bnds.f = [1 500];
bnds.k11 = [1e-5 200];
bnds.k22 = [0 1];
bnds.k21 = [1 200];
bnds.kF = [1 1e5];
bnds.J1 = [1e-3 200];
bnds.J2 = [1 1e3];
bnds.kon = [1 200];
bnds.kse = [1e-3 1];
bnds.kse0 = [1e-4 1];
bnds.koop = [1 200];
bnds.n = [.1 10];
bnds.kappa = [.1 10];
bnds.vmax = [1 200];
bnds.ps2 = [-1 2];
bnds.k = [1 5000];
bnds.b = [1 5000];
bnds.dLcrit = [1 4];

% parameters to be fitted
if sum(mcode == [1 1 1]) == 3
    optparms = {'f', 'k11', 'k22', 'k21', 'kon', 'kse', 'kse0', 'kF', 'koop'};
%     optparms = {'f', 'k11', 'k22', 'k21', 'kon', 'kse', 'kse0', 'kF'};
elseif sum(mcode == [1 1 3]) == 3
    optparms = {'f', 'k11', 'k22', 'k21', 'n','kappa', 'kse', 'kse0'};
elseif sum(mcode == [2 1 1]) == 3
    optparms = {'n','kappa', 'kse', 'kse0', 'vmax'};
elseif sum(mcode == [1 2 1]) == 3
    optparms = {'f', 'k11', 'k22', 'k21', 'kon', 'kse', 'kse0', 'kF', 'koop','ps2', 'dLcrit'};
end

%% specify data
load('active_trials.mat', 'Fm')

for iF = iFs
    
    %% load data
    cd(['C:\Users\',username,'\OneDrive - KU Leuven\9. Short-range stiffness\matlab\data'])
    load([fibers{iF},'_cor_new.mat'],'data')
   
    Ks = find(Fm(:,iF) > .07); % only consider active trials
    Data = prep_data_v2(data,n, m,Ks,tiso);
    [tis, Cas, Lis, vis, ts] = create_input(tiso, Data.dTt, Data.dTc, Data.ISI, Data.Ca(Ks), 500);
    
    Lis(Lis<0) = 0;
    
    % interpolate force
    Fis = interp1(Data.t, Data.F, tis);
    
    if ishandle(1), close(1); end
    figure(1)
    subplot(411)
    plot(Data.t, Data.C,'k.', tis, Cas, 'b');
    
    subplot(412)
    plot(Data.t, Data.v,'k.', tis, vis, 'b');
    
    subplot(413)
    plot(Data.t, Data.L,'k.', tis, Lis,'b');
    
    subplot(414)
    plot(Data.t, Data.F,'k.');
    
    for j = 1:4
        subplot(4,1,j)
        box off
        hold on
    end
    
    %% get parameters
    vs = {'\'};
  
    if sum(mcode == [1 2 1]) == 3 % use parms from coop model
        mmcode = [1 1 1];
    else
        mmcode = mcode; % use own parms
    end
    
    cd([githubfolder, '\muscle-thixotropy\new_model\get_variable'])
    [output_mainfolder, filename, opt_type, ~] = get_folder_and_model(mmcode);
    
    disp(filename)
%     output_folder = [opt_type,'\normalized\with_PE_optimized\2_trials'];
    
    % get other parameters from other fiber    
%     output_dir = [output_mainfolder{1}, '\', filename,vs{1}, output_folder];
%     cd(output_dir)
%     load([filename,'_F', num2str(iF),'_best.mat'],'parms','exitflag','fopt','C0','Cbounds','model','P0','P')
    
    foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iF}];
    cd(foldername)
    load(['parms_', filename, '_v2.mat'], 'newparms')

    parms = newparms;
    % convert to parms
%     C = p_to_c(P, Cbounds);
%     parms = C_to_parms(C, parms, parms.optvars);
%     parms = calc_dependent_parms(parms);
    
    % add some parameters that didn't exist before
%     parms.ps2 = 0;
    parms.act = 1;
    parms.Noverlap = 1;
    parms.approx = 1;
    parms.K = 100;
    parms.vF_func = @(vcerel,parms)parms.e(1)*log((parms.e(2)*vcerel./parms.vmax+parms.e(3))+sqrt((parms.e(2)*vcerel./parms.vmax+parms.e(3)).^2+1))+parms.e(4);
    
    if ~isfield(parms, 'gamma')
        parms.gamma = .5*parms.s / parms.h; % length scaling
    end
    
    % store old parameters
    oparms = parms;
    
    % update powerstroke
    h = 10e-9; % powerstroke size
    s = 2.6e-6; % sarcomere length
    parms.gamma = .5*s / h; % length scaling
    parms.kpe = oparms.gamma / parms.gamma * oparms.kpe;
    
    if sum(mcode == [1 1 1]) == 3
        parms.J1 = 6.17;
        parms.koop = 5.7;
        parms.JF = 1e3;
        parms.J2 = 200;
        
    elseif sum(mcode == [2 1 1]) == 3 % Hill-type
        parms.f = 0;
        parms.k = 0;
        parms.dLcrit = 0;
        
        parms.act_max = (1 - parms.Fpe0) / parms.Fscale;

    elseif sum(mcode == [1 2 1]) == 3 
%         parms.J1 = 6.17;
%         parms.koop = 5.7;
%         parms.JF = 1e3;
        
        parms.k = 3000;
        parms.b = 1000;
        parms.dLcrit = 2.5;
    end
    
    %% get initial guess
    parms.ti = tis;
    parms.vts = vis;
    parms.Cas = Cas;
    parms.Lts = Lis * parms.gamma;
    
    IG = get_initial_guess(tis, Cas, vis, parms.Lts, parms);
    
    oFi = IG.Fi * parms.Fscale + parms.Fpe_func(parms.Lts, parms);
    
    subplot(414); hold on
    plot(tis, oFi,'b'); hold on
    
    %% intervals of interest
    id1 = nan(length(Ks), length(parms.ti), length(m));
    id2 = nan(length(Ks), length(parms.ti), length(m));
    tmax = [0 tiso * length(Ks)];
    
    for i = 1:length(m)
        for k = 1:length(Ks)
            id1(k,:,i) = (parms.ti > (tmax(i) + ts(k) + tiso - 4*Data.dTt - 2*Data.dTc(i) - Data.ISI(i))) & (parms.ti < (tmax(i) + ts(k)+tiso - Data.dTt));
            
            id2(k,:,i) = (parms.ti > (tmax(i) + ts(k) + tiso - 5*Data.dTt - 2*Data.dTc(i) - Data.ISI(i))) & (parms.ti < (tmax(i) + ts(k) + tiso - 4*Data.dTt - 2*Data.dTc(i) - Data.ISI(i)));
        end
    end
    
    % sum across conditions and pCas
    idF = find(sum(sum(id1,3, 'omitnan'),1) & isfinite(Fis));
    idC = find(sum(sum(id2,3, 'omitnan'),1));

    subplot(414); hold on
    plot(tis(idF), oFi(idF),'m*', 'markersize', 1); hold on
    plot(tis(idC), oFi(idC),'g.'); hold on
    
    %% select data for fitting
    Xdata.t = tis;
    Xdata.F = Fis;
    Xdata.v = vis;
    Xdata.L = Lis * parms.gamma;
    Xdata.Cas = Cas;
    Xdata.idF = idF;
    Xdata.idC = idC;
    
    %% do fitting    
    % initialise opti structure
    opti = casadi.Opti();
    
    % define weigth vector
    w1 = 100;    % weight for fitting force-velocity
    w2 = 100;   % weight for fitting short-range stiffness
    w3 = 1000; 	% weight for regularization
    w = [w1 w2 w3];
    
    % specify biophysical parameters to be fitted
    parms.kF = parms.J1 * parms.JF;
    fparms = parms;
    
    figure(100)
    [newparms, out] = fit_model_parameters_v2(opti, optparms, w, Xdata, fparms, IG, bnds);
    set(gcf,'units','normalized','position',[.2 .2 .4 .6])
    
    if newparms.J1 > 0
        newparms.JF = newparms.kF / newparms.J1;
    end

    %% test with fitted paramers   
    oLiss = Lis * oparms.gamma;
    nLiss = Lis * newparms.gamma;
    
    oparms.ti = tis;
    oparms.vts = vis;
    oparms.Cas = Cas;
    oparms.Lts = oLiss;
    
    newparms.ti = tis;
    newparms.vts = vis;
    newparms.Cas = Cas;
    newparms.Lts = nLiss;
    
    odeopt = odeset('maxstep', 3e-3);
    
    if newparms.f > 0
        x0 = 1e-3 * ones(7,1);
        
        if newparms.J1 == 0
            x0(6) = 1;
        end
        
        xp0 = zeros(size(x0));
        
        osol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, oparms), [0 max(tis)], x0, xp0, odeopt);
        nsol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, newparms), [0 max(tis)], x0, xp0, odeopt);
        nF = (nsol.y(1,:) + nsol.y(2,:)) * parms.Fscale;
        oF = (osol.y(1,:) + osol.y(2,:)) * oparms.Fscale;
        
        nt = nsol.x;
        ot = osol.x;
        
        nFi = interp1(nt, nF, tis) + parms.Fpe_func(nLiss, newparms);
        oFi = interp1(ot, oF, tis) + parms.Fpe_func(oLiss, oparms);
        
    else
        
        x0 = 0;
        xp0 = zeros(size(x0));
        
        % simulate
        osol = ode15i(@(t,y,yp) hill_type_implicit_v2(t,y,yp, oparms), [0 max(tis)], x0, xp0, odeopt);
        nsol = ode15i(@(t,y,yp) hill_type_implicit_v2(t,y,yp, newparms), [0 max(tis)], x0, xp0, odeopt);
        
        % elastic elements
        oL = oLiss - interp1(osol.x, osol.y(1,:), tis);
        nL = nLiss - interp1(nsol.x, nsol.y(1,:), tis);
        
        oFi = parms.Fse_func(oL, oparms) * parms.Fscale + parms.Fpe_func(oLiss, oparms);
        nFi = parms.Fse_func(nL, newparms) * parms.Fscale + parms.Fpe_func(nLiss, newparms);
        
    end
    
    %% compare forces    
    if ishandle(2), close(2); end
    
    figure(2)
    plot(tis(Xdata.idF), Fis(idF),'k.'); hold on
    plot(out.t, out.F, 'linewidth', 2)
    plot(tis, nFi,'--', 'linewidth', 2); 
    plot(tis, oFi,'--', 'linewidth', 2);
    
    box off
    legend('Data', 'Fitted', 'New', 'Old', 'location', 'best')
    legend boxoff
    
    %% parameter values
    figure(3)
    nexttile
    bar(categorical(optparms), out.s)
    
    %% SRS state    
%     figure(4)
%     nexttile
%     plot(osol.x, osol.y(end,:)); hold on
%     plot(nsol.x, nsol.y(end,:))
    
    %% estimate SRS
    nFii = interp1(tis, nFi, Data.t);
    oFii = interp1(tis, oFi, Data.t);
       
    if ishandle(5), close(5); end
    figure(5)
    subplot(131);
    plot(Data.L, Data.F,'.'); hold on
    
    subplot(132);
    plot(Data.L, oFii,'.'); hold on
    
    subplot(133)
    plot(Data.L, nFii,'.'); hold on
    
    % pre-allocate
    ns = nan(length(Ks), 2);
    os = nan(length(Ks), 2);
    ds = nan(length(Ks), 2);
    F0 = nan(length(Ks), 1);
    
    [id0,id1,id2] = get_indices(Data.t, tiso, ts, Data.dTt, Data.dTc(1), Data.ISI(1), Data.Ca(Ks));
    
    for i = 1:length(Ks)
        
        np1 = polyfit(Data.L(id1(i,:)), nFii(id1(i,:)), 1);
        np2 = polyfit(Data.L(id2(i,:)), nFii(id2(i,:)), 1);
        
        op1 = polyfit(Data.L(id1(i,:)), oFii(id1(i,:)), 1);
        op2 = polyfit(Data.L(id2(i,:)), oFii(id2(i,:)), 1);
        
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
        plot(Data.L(id1(i,:)), oFii(id1(i,:)), 'r.')
        plot(Data.L(id2(i,:)), oFii(id2(i,:)), 'y.')
        
        subplot(133)
        plot(Data.L(id1(i,:)), nFii(id1(i,:)), 'r.')
        plot(Data.L(id2(i,:)), nFii(id2(i,:)), 'y.')
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
    figure(6)
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

    %% save
    if save_results

       [output_mainfolder, filename, opt_type, ~] = get_folder_and_model(mcode);

        foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iF}];
        if ~isfolder(foldername)
            mkdir(foldername)
        end

        cd(foldername)
        save(['parms_', filename, '_v2.mat'], 'newparms', 'optparms', 'out', 'bnds')
    end
    
end


