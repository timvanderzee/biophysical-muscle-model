clear all; close all; clc
[username, githubfolder] = get_paths();
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

% model to be fitted
mcode = [1 1 1];
version = '_v3'; 
% version = '';

% settings
N = 500;
% N = 1500;
save_results = 0;
visualize = 1;
    
iFs = 11 %; [2,3,5,6,7,8,11];
% n = [3 1]; % ISI number
% m = [7 1]; % AMP number
n = [3]; % ISI number
m = [7]; % AMP number
tiso = 3; % isometric time (s)

% bounds
bnds.f = [20 500];
bnds.k11 = [1e-5 200];
bnds.k22 = [0 1];
bnds.k21 = [1 200];
bnds.kF = [1 1e4];
bnds.J1 = [1e-3 200];
bnds.J2 = [1 1e3];
bnds.kon = [5 200];
bnds.kse = [1e-3 1];
% bnds.kse0 = [1e-4 1];
bnds.kse0 = [1e-2 1];
bnds.koop = [1 200];
bnds.n = [.1 10];
bnds.kappa = [.1 10];
bnds.vmax = [1 200];
bnds.ps2 = [-1 2];
bnds.k = [1 5000];
bnds.b = [1 5000];
bnds.dLcrit = [1 10];
bnds.Lce0 = [-10 0];

bnds.kpe = [1e-4 1e-1];
bnds.Fpe0 = [1e-5 1e-1];

% parameters to be fitted
if sum(mcode == [1 1 1]) == 3
    optparms = {'f', 'k11', 'k22', 'k21', 'kon', 'kse', 'kse0', 'kF', 'koop', 'Lce0', 'kpe'};
%     optparms = {'kpe', 'Fpe0'};
%     optparms = {'f', 'k11', 'k22', 'k21', 'kon', 'kse', 'kse0', 'kF'};
elseif sum(mcode == [1 1 2]) == 3
    optparms = {'f', 'k11', 'k22', 'k21', 'kon', 'kse', 'kse0', 'koop', 'kpe', 'Fpe0'};
elseif sum(mcode == [1 1 3]) == 3
    optparms = {'f', 'k11', 'k22', 'k21', 'n','kappa', 'kse', 'kse0'};
elseif sum(mcode == [2 1 1]) == 3
    optparms = {'n','kappa', 'kse', 'kse0', 'vmax'};
elseif sum(mcode == [1 2 1]) == 3
    optparms = {'f', 'k11', 'k22', 'k21', 'kon', 'kse', 'kse0', 'kF', 'koop', 'dLcrit'};
end

%% specify data
load('active_trials.mat', 'Fm')

for iF = iFs
    
    %% load data
    cd(['C:\Users\',username,'\OneDrive - KU Leuven\9. Short-range stiffness\matlab\data'])
    load([fibers{iF},'_cor_new.mat'],'data')
   
    Ks = find(Fm(:,iF) > 0); % only consider active trials
    Data = prep_data_v2(data,n, m,Ks,tiso);
    [tis, Cas, Lis, vis, ts] = create_input(tiso, Data.dTt, Data.dTc, Data.ISI, Data.Ca(Ks), N);
    
    Lis(Lis<0) = 0;
    
    % interpolate force
    Fis = interp1(Data.t, Data.F, tis);
    
    if visualize
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
    end
    
    %% get parameters
    vs = {'\'};
  
    if sum(mcode == [1 2 1]) == 3 % use parms from coop model
        mmcode = [1 1 1];
    elseif sum(mcode == [1 1 2]) == 3
        mmcode = [1 1 3];
    else
        mmcode = mcode; % use own parms
    end
    
    cd([githubfolder, '\muscle-thixotropy\new_model\get_variable'])
    [output_mainfolder, filename, opt_type, ~] = get_folder_and_model(mmcode);
    
    disp(filename)
    
    foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iF}];
    cd(foldername)
    load(['parms_', filename, version, '.mat'], 'newparms')

    parms = newparms;
    
    % add some parameters that didn't exist before
    parms.PE_isw_SE = 1;
    parms.Fse_func = @(dlse, parms) parms.kse0*(exp(parms.kse*dlse)-1);
    parms.Fpe0 = parms.Fpe0/2;
    parms.Lce0 = -1;
    parms.Fpe_func = @(L, parms) parms.kpe*(L-parms.lmtc0).*(L>parms.lmtc0)+parms.Fpe0;
    parms.kpe_func = @(Lce, parms) parms.kpe;

    parms.kF = min(parms.kF, .9*bnds.kF(2));
    parms.kse0 = max(parms.kse0, 1.1*bnds.kse0(1));
    
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
%         parms.J1 = 6.17;
%         parms.koop = 5.7;
%         parms.JF = 1e3;
%         parms.J2 = 200;
        
    elseif sum(mcode == [1 1 2]) == 3 % thin coop only
        parms.actfunc = @(Ca,parms) Ca;
        parms.koop = 5;
        parms.koff = 80;
        parms.kon  = 33;
        
    elseif sum(mcode == [2 1 1]) == 3 % Hill-type
        parms.f = 0;
        parms.k = 0;
        parms.dLcrit = 0;
        
        parms.act_max = (1 - parms.Fpe0) / parms.Fscale;

    elseif sum(mcode == [1 2 1]) == 3 
        parms.J1 = 6.17;
        parms.koop = 5.7;
        parms.JF = 1e3;
        parms.J2 = 200;
        
        parms.k = 1000;
        parms.b = 3000;
        parms.dLcrit = 2;
        
        parms.ps2 = 0;
    end
    
    %% get initial guess
    parms.ti = tis;
    parms.vts = vis;
    parms.Cas = Cas;
    parms.Lts = Lis * parms.gamma;
    
    IG = get_initial_guess(tis, Cas, vis, parms.Lts, parms);
    
    oFi = (IG.Fsei) * parms.Fscale;
    
    if visualize
        figure(1)
        subplot(414); hold on
        plot(tis, oFi,'b'); hold on
    end
    
    Lcost = @(dlse, Lis, parms) ((parms.kse0*(exp(parms.kse*dlse)-1)) - (parms.kpe * (Lis*parms.gamma - dlse - parms.Lce0))).^2;

    for i = 1:length(Lis)
        dlse(i) = fminsearch(@(L) Lcost(L, Lis(i), parms), 0);
    end

    dLce = (Lis * parms.gamma - dlse) - parms.Lce0;

    Fse = parms.Fse_func(dlse, parms) * parms.Fscale;
    Fpe = parms.kpe * dLce  * parms.Fscale;

    figure(1)
    subplot(414)
    plot(tis, Fse, '-', tis, Fpe, '--')

    % initial guess
    IG.F9 = Fse;

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

    if visualize
        subplot(414); hold on
        plot(tis(idF), oFi(idF),'m*', 'markersize', 1); hold on
        plot(tis(idC), oFi(idC),'g.'); hold on
    end
    
    % active indices (assume last trial is passive)
    idA = 1:length(tis) - N;
    
    %% select data for fitting
    Xdata.t = tis;
    Xdata.F = Fis;
    Xdata.v = vis;
    Xdata.L = Lis * parms.gamma;
    Xdata.Cas = Cas;
    Xdata.idF = idF(idF<idA(end));
    Xdata.idC = idC(idC<idA(end));
    
    % active
    Xdata.idA = idA;
    
    % passive
    Xdata.idP = (idA(end)+1):length(tis);
    Xdata.idFP = idF(idF>idA(end));
    
    %% do fitting        
    % define weigth vector
    w1 = 100;    % weight for fitting force-velocity
    w2 = 1000;   % weight for fitting short-range stiffness
    w3 = 1000; 	% weight for regularization
    w = [w1 w2 w3];
    
    % specify biophysical parameters to be fitted
%     parms.kF = parms.J1 * parms.JF;
%     parms.Lce0 = -20;
    fparms = parms;

%     figure(100)
    [newparms, out, opti] = fit_model_parameters_v2(optparms, w, Xdata, fparms, IG, bnds);
    set(gcf,'units','normalized','position',[.2 .2 .4 .6])
    
    if newparms.J1 > 0
        newparms.JF = newparms.kF / newparms.J1;
    end
    
     
    %% test with fitted paramers     
    close all
    N = 1000;
    [tis, Cas, Lis, vis, ts] = create_input(tiso, Data.dTt, Data.dTc, Data.ISI, Data.Ca(Ks), N);
     
    nLiss = Lis * newparms.gamma;

    newparms.ti = tis;
    newparms.vts = vis;
    newparms.Cas = Cas;
    newparms.Lts = nLiss;
    
    odeopt = odeset('maxstep', 1e-2);
    
    % initial conditions
    Q0 = .1;
    p0 = -.5;
    q0 = .1;
    [Q00, Q20, lce0, Q10, Fse0] = find_steady_state(Q0, p0, q0, newparms, 'regular');
    x0 = [Q00 Q20 Fse0 0 0 0];

    nsol = ode15s(@(t,y) fiber_dynamics_explicit_length_v2(t,y, newparms), [0 max(tis)], x0, odeopt);
    nt = nsol.x;

    % calc force
    nF = nsol.y(3,:);
    ndlse = parms.Lse_func(nF, newparms);
    L = newparms.Lts - interp1(nt, ndlse, newparms.ti);
    nFi = interp1(nt, nF, tis) * parms.Fscale; % + nFp;
    
    if ishandle(2), close(2); end
    figure(2)
    plot(Data.t, Data.F,'k.'); hold on
    plot(out.t, out.F)
    plot(tis, nFi,'--');

    legend('Data', 'Fit', 'Sim')

    %% estimate SRS
    % recreate the input, but with a higher sampling rate

    nFii = interp1(tis, nFi, Data.t);

    Lti = interp1(tis, newparms.Lts, Data.t);
    
    % pre-allocate
    ns = nan(length(Ks), 2);
    ds = nan(length(Ks), 2);
    F0 = nan(length(Ks), 1);
    
    [id0,id1,id2] = get_indices(Data.t, tiso, ts, Data.dTt, Data.dTc(1), Data.ISI(1), Data.Ca(Ks));
    
    for i = 1:length(Ks)
        
        np1 = polyfit(Lti(id1(i,:)), nFii(id1(i,:)), 1);
        np2 = polyfit(Lti(id2(i,:)), nFii(id2(i,:)), 1);
        
        dp1 = polyfit(Data.L(id1(i,:)), Data.F(id1(i,:)), 1);
        dp2 = polyfit(Data.L(id2(i,:)), Data.F(id2(i,:)), 1);
        
        ds(i,:) = [dp1(1) dp2(1)];
        ns(i,:) = [np1(1) np2(1)];
        
        F0(i) = mean(Data.F(id0(i,:)));
        
    end

    
    figure(6)
    plot(F0, ds(:,1)./ds(:,2),'o-'); hold on
    plot(F0, ns(:,1)./ns(:,2),'o-')
    box off
    xlabel('Isometric force (F_0)')
    ylabel('Relative stiffness')
    legend('Data', 'Model', 'location', 'best')
    legend boxoff
    xlim([0 1.05])
    
    return
     
    %% load 
    [output_mainfolder, filename, opt_type, ~] = get_folder_and_model(mmcode);

    foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iF}];
    if ~isfolder(foldername)
        mkdir(foldername)
    end

    cd(foldername)

    oldparms = load(['parms_', filename, version, '.mat'], 'newparms', 'optparms', 'out', 'bnds');
    oparms = oldparms.newparms;
    
    %% test with fitted paramers 
    N = 5e3;
    [tis, Cas, Lis, vis, ts] = create_input(tiso, Data.dTt, Data.dTc, Data.ISI, Data.Ca, N);
      
    oparms.approx = 0;
    newparms.approx = 0;
    
    oLiss = Lis * oparms.gamma;
    nLiss = Lis * newparms.gamma;
    
%     oparms.Fpe_func = @(Lce, parms) parms.kpe * log(1+exp(Lce*parms.K))/parms.K + parms.Fpe0;
%     oparms.kpe_func = @(Lce, parms) parms.kpe .* (1 - 1./(exp(parms.K*Lce)+1));
    
    oparms.ti = tis;
    oparms.vts = vis;
    oparms.Cas = Cas;
    oparms.Lts = oLiss;
    
    newparms.ti = tis;
    newparms.vts = vis;
    newparms.Cas = Cas;
    newparms.Lts = nLiss;
    
    newparms.kpe = 0;
    newparms.Fpe0 = 0;
    
    odeopt = odeset('maxstep', 1e-2);
    
    % initial conditions
    Q00 = 0;
    p0 = 0;
    q0 = .1;
%     [Q00, Q20, lce0, Q10] = find_steady_state(Q0, p0, q0, parms);

    if newparms.f > 0
%         x0 = 1e-3 * ones(7,1);
        
        if newparms.J1 == 0
            x0(6) = 1;
        end
        
        xp0 = zeros(size(x0));
        
        [Q00, Q20, lce0, Q10] = find_steady_state(Q00, p0, q0, oparms, 'adjusted');
        x0 = [Q00 Q10 Q20 lce0 0 0 0];
        osol = ode15s(@(t,y) fiber_dynamics_explicit_no_tendon(t,y, oparms), [0 max(tis)], x0, odeopt);
        
        Q0 = osol.y(1,:);
        Q1 = osol.y(2,:);
        oF = (Q0 + Q1) * parms.Fscale + oparms.Fpe_func(interp1(oparms.ti, oparms.Lts, osol.x), oparms);
        ot = osol.x;
        oFi = interp1(ot, oF, tis);
       

        [Q00, Q20, lce0, Q10] = find_steady_state(Q00, p0, q0, newparms, 'regular');
        x0 = [Q00 Q20 lce0 0 0 0];
        nsol = ode15s(@(t,y) fiber_dynamics_explicit_length(t,y, newparms), [0 max(tis)], x0, odeopt);
        nt = nsol.x;
        
        % calc force
        nL = nsol.y(end-3,:);
        ndlse = interp1(newparms.ti, newparms.Lts, nsol.x) - nL;
        nF = newparms.Fse_func(ndlse, newparms) * parms.Fscale;
        
        nFp = newparms.Fpe_func(nL, newparms);
        
        nFi = interp1(nt, nF, tis);
 
        
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
    Data = prep_data_v2(data,n, m,1:7,tiso);
        
    if visualize
        if ishandle(2), close(2); end

        figure(2)
        plot(Data.t, Data.F,'k.'); hold on
%         plot(out.t, out.F, 'linewidth', 2)
        plot(tis, nFi,'-', 'linewidth', 2); 
        plot(tis, oFi,'-', 'linewidth', 2);

        box off
        legend('Data', 'New', 'Old', 'location', 'best')
        legend boxoff
    end
    
    %% parameter values
    for i = 1:length(optparms)
        lb(i) = bnds.(optparms{i})(1);
        ub(i) = bnds.(optparms{i})(2);
        nv(i) = (oldparms.newparms.(optparms{i}) - lb(i)) / (ub(i)-lb(i));
    end
    
    figure(3)
    nexttile
   
    bar(categorical(optparms), [out.s; nv])
    
    %% SRS state    
%     figure(4)
%     nexttile
%     plot(osol.x, osol.y(end,:)); hold on
%     plot(nsol.x, nsol.y(end,:))
    

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
        save(['parms_', filename, '_v3.mat'], 'newparms', 'optparms', 'out', 'bnds')
    end
    
end

