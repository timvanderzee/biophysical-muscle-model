clear all; close all; clc
[username, githubfolder] = get_paths();
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

% model to be fitted
mcode = [1 1 1];
old_version = '_v3';
new_version = '_v4';

% settings
N = 500;
% N = 1500;
save_results = 0;
visualize = 1;
    
iFs = 2 %; [2,3,5,6,7,8,11];
n = [3 1]; % ISI number
m = [7 1]; % AMP number
% n = [1]; % ISI number
% m = [1]; % AMP number
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
    optparms = {'f', 'k11', 'k22', 'k21', 'kon', 'kse', 'kse0', 'koop','Lce0', 'kpe'};
elseif sum(mcode == [1 1 3]) == 3
    optparms = {'f', 'k11', 'k22', 'k21', 'n','kappa', 'kse', 'kse0','Lce0', 'kpe'};
elseif sum(mcode == [2 1 1]) == 3
    optparms = {'n','kappa', 'kse', 'kse0', 'vmax'};
elseif sum(mcode == [1 2 1]) == 3
    optparms = {'f', 'k11', 'k22', 'k21', 'kon', 'kse', 'kse0', 'kF', 'koop', 'dLcrit', 'Lce0', 'kpe'};
end

%% specify data
load('active_trials.mat', 'Fm')

%%

for iF = iFs
    
    %% load data
    cd(['C:\Users\',username,'\OneDrive - KU Leuven\9. Short-range stiffness\matlab\data'])
    load([fibers{iF},'_cor_new.mat'],'data')
   
    Ks = find(Fm(:,iF) > .05); % only consider active trials
    AData = prep_data_v2(data,n,m,Ks,tiso);
    [Atis, ACas, ALis, Avis, Ats] = create_input(tiso, AData.dTt, AData.dTc, AData.ISI, AData.Ca(Ks), N);
    
    % passive
    PData = prep_data_v2(data,3,7,length(data.pCas),tiso);
    [Ptis, PCas, PLis, Pvis, Pts] = create_input(tiso, PData.dTt, PData.dTc, PData.ISI, PData.Ca(length(data.pCas)), N);
    
    % combine
    Data.t = [AData.t; PData.t(2:end) + AData.t(end)];
    Data.F = [AData.F; PData.F(2:end)];
    Data.L = [AData.L; PData.L(2:end)];
    Data.v = [AData.v; PData.v(2:end)];
    Data.C = [AData.C; PData.C(2:end)];
    
    % combine model input
    tis = [Atis Ptis(2:end) + Atis(end)];
    Cas = [ACas PCas(2:end)];
    Lis = [ALis PLis(2:end)];
    vis = [Avis Pvis(2:end)];
    ts = [Ats Pts + Ats(end)];

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
        mmcode = [1 1 1];
    else
        mmcode = mcode; % use own parms
    end
    
    cd([githubfolder, '\muscle-thixotropy\new_model\get_variable'])
    [output_mainfolder, filename, opt_type, ~] = get_folder_and_model(mmcode);
    
    disp(filename)
    
    foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iF}];
    cd(foldername)
    load(['parms_', filename, old_version, '.mat'], 'newparms')

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
        parms.J1 = 0;
        parms.J2 = 0;
        parms.JF = 0;
        parms.kF = 0;
        
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
    
    %% intervals of interest
    parms.ti = tis;
    parms.vts = vis;
    parms.Cas = Cas;
    parms.Lts = Lis * parms.gamma;
    
    id1 = nan(length(Ks), length(parms.ti), length(m));
    id2 = nan(length(Ks), length(parms.ti), length(m));
    tmax = [0 tiso * length(Ks) tiso * (length(Ks)*2)];
    
    for i = 1:length(m)
        for k = 1:length(Ks)
            id1(k,:,i) = (parms.ti > (tmax(i) + ts(k) + tiso - 4*AData.dTt - 2*AData.dTc(i) - AData.ISI(i))) & (parms.ti < (tmax(i) + ts(k)+tiso - AData.dTt));
            id2(k,:,i) = (parms.ti > (tmax(i) + ts(k) + tiso - 5*AData.dTt - 2*AData.dTc(i) - AData.ISI(i))) & (parms.ti < (tmax(i) + ts(k) + tiso - 4*AData.dTt - 2*AData.dTc(i) - AData.ISI(i)));
        end
    end
    
    % sum across conditions and pCas
    idF = find(sum(sum(id1,3, 'omitnan'),1) & isfinite(Fis));
    idC = find(sum(sum(id2,3, 'omitnan'),1));
    
    % active indices (assume last trial is passive)
    idA = 1:length(tis) - N;
    
    % passive
    idP = (idA(end)+1):length(tis);
    idFP = find((parms.ti > (tmax(3) + tiso - 4*PData.dTt - 2*PData.dTc - PData.ISI)) & (parms.ti < (tmax(3) + tiso - PData.dTt)) & isfinite(Fis));
      
    
    %% get initial guess
    IG = get_initial_guess(tis, Cas, vis, parms.Lts, parms);
    
    oFi = (IG.Fsei) * parms.Fscale;
    
    if visualize
        figure(1)
        subplot(414); hold on
        plot(tis, oFi,'b'); hold on
        
        plot(tis(idF), oFi(idF),'m*', 'markersize', 1); hold on
        plot(tis(idC), oFi(idC),'g.'); hold on
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
    Xdata.idP = idP;
    Xdata.idFP = idFP;
    
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
    [newparms, out, opti] = fit_model_parameters_v2(optparms, w, Xdata, fparms, IG, bnds, mcode);
    set(gcf,'units','normalized','position',[.2 .2 .4 .6])
    
    if newparms.J1 > 0
        newparms.JF = newparms.kF / newparms.J1;
    end
    
    %% sim output
    Ks = find(Fm(:,iF) > .05);
    Data = prep_data_v2(data,n,m,Ks,tiso);

    if ishandle(2), close(2); end
    figure(2)
    subplot(121);
    plot(Data.t, Data.F,'k.'); hold on
    plot(out.t, out.F) 

    for i = 1:length(optparms)
        lb(i) = bnds.(optparms{i})(1);
        ub(i) = bnds.(optparms{i})(2);
    end
    
    subplot(122);
    bar(categorical(optparms), out.s)
     
    %% test with fitted paramers     
    Ks = find(Fm(:,iF) > 0);
    Data = prep_data_v2(data,n,m,Ks,tiso);

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
    p0 = -1;
    q0 = .1;
    [Q00, Q20, lce0, Q10, Fse0] = find_steady_state(Q0, p0, q0, newparms, 'regular');
    x0 = [Q00 Q20 Fse0 0 0 0];
    x0 = zeros(1,6);
    x0(5) = 1;
    
    nsol = ode15s(@(t,y) fiber_dynamics_explicit_length_v2(t,y, newparms), [0 max(tis)], x0, odeopt);
    nt = nsol.x;

    % calc force
    nF = nsol.y(3,:);
    ndlse = parms.Lse_func(nF, newparms);
    L = newparms.Lts - interp1(nt, ndlse, newparms.ti);
    nFi = interp1(nt, nF, tis) * parms.Fscale; % + nFp;
    
    if ishandle(3), close(3); end
    figure(3)
    plot(Data.t, Data.F,'k.'); hold on
    plot(tis, nFi,'-');

    legend('Data', 'Sim')

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

    
    if ishandle(6), close(6); end; figure(6)
    plot(F0, ds(:,1)./ds(:,2),'o-'); hold on
    plot(F0, ns(:,1)./ns(:,2),'o-')
    box off
    xlabel('Isometric force (F_0)')
    ylabel('Relative stiffness')
    legend('Data', 'Model', 'location', 'best')
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
        save(['parms_', filename, new_version,'.mat'], 'newparms', 'optparms', 'out', 'bnds')
    end
    
end

