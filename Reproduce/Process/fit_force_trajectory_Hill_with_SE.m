clear all; close all; clc
[username, githubfolder] = get_paths();
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

% model to be fitted
mcode = [2 2 1];
old_version = '';
new_version = '';

% settings
N = 500;
% N = 1500;
save_results = 0;
visualize = 1;
    
iFs = 6 %; [2,3,5,6,7,8,11];
n = [3 1]; % ISI number
m = [7 1]; % AMP number
% n = [1]; % ISI number
% m = [1]; % AMP number
tiso = 3; % isometric time (s)


%% specify data
load('active_trials.mat', 'Fm')


for iF = iFs
    
    %% load data
    cd(['C:\Users\',username,'\OneDrive - KU Leuven\9. Short-range stiffness\matlab\data'])
    load([fibers{iF},'_cor_new.mat'],'data')
   
    Ks = find(Fm(:,iF) > 0); % only consider active trials
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
  
    if isequal(mcode, [2 2 1])
        mmcode = [2 1 1]; % use own parms
    end
    
    cd([githubfolder, '\muscle-thixotropy\new_model\get_variable'])
    [output_mainfolder, filename, opt_type, ~] = get_folder_and_model(mmcode);
    
    disp(filename)
    
    foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iF}];
    cd(foldername)
    load(['parms_', filename, old_version, '.mat'], 'newparms')

    parms = newparms;
 
    
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

    
    %% get initial guess       
    parms.act_max = .5;
    
    a = parms.actfunc(Cas, parms);
    a(a<parms.amin) = parms.amin;
    
    parms.Lce0 = 0;

    IG.Fcei = parms.vF_func(vis * parms.gamma, parms) .* a;
    IG.Fpei = parms.kpe * (Lis - parms.Lce0) * parms.gamma .* (Lis > parms.Lce0);

    IG.Fsei = IG.Fcei + IG.Fpei;
    
    %% fit
    import casadi.*

    % initialise opti structure
    opti = casadi.Opti();
    
    % parameters
    kappa   = opti.variable(1);
    act_max = opti.variable(1);
    vmax    = opti.variable(1);
    n       = opti.variable(1);
    kpe     = opti.variable(1);
    Lce0    = opti.variable(1);
    kse0     = opti.variable(1);
    kse      = opti.variable(1);
    
    % initial values
    opti.set_initial(kappa, parms.kappa);
    opti.set_initial(act_max, parms.act_max);
    opti.set_initial(vmax, parms.vmax);
    opti.set_initial(n, parms.n);
    opti.set_initial(kpe, parms.kpe);
    opti.set_initial(Lce0, parms.Lce0);
    opti.set_initial(kse0, parms.kse0);
    opti.set_initial(kse, parms.kse);
    
    % bounds
    opti.subject_to(kappa > 0);
    opti.subject_to(act_max > 0);
    opti.subject_to(vmax > 0);
    opti.subject_to(n > 0);
    opti.subject_to(kpe > 0);
    opti.subject_to(Lce0 < 0);
    opti.subject_to(kse0 > 0);
    opti.subject_to(kse > 0);
    
    % variables
    Lts = Lis * parms.gamma;
    vts = vis * parms.gamma;
    Fts = Fis;
    
    N = length(Lis);
    F = opti.variable(1, N);
    
    dlse = log(F);
    
    L = Lts - dlse;
    Fpe = F(L);
    
    
    Fce = F - Fpe;
    v = f(Fce)
    
    Fdot = f(v);
    
    
    % PE
    dLce = L - Lce0;
    K = 1;
    Fpe = kpe * log(1+exp(dLce*K))/K;
    
    % CE
    a = act_max * Cas.^n ./ (kappa^n + Cas.^n);
    Fce = a .* (parms.e(1)*log((parms.e(2)*vts./vmax+parms.e(3))+sqrt((parms.e(2)*vts./vmax+parms.e(3)).^2+1))+parms.e(4));
    
    % SE
    dlse = Lts - L;
    Fse = kse0 * (exp(kse*dlse)-1);
    
    % Total force
    Frel = (Fce + Fpe) * parms.Fscale;
    
    opti.subject_to(Frel(1) == 1);

    Fcost = (Frel(idF) - Fts(idF)).^2;

    % cost function
    J = 0;
    J = J + sum(Fcost);

    % optimize
    opti.minimize(J);

    % options for IPOPT
    options.ipopt.linear_solver     = 'mumps';
    options.ipopt.mu_strategy       = 'adaptive';
    options.detect_simple_bounds    = true;
    options.ipopt.max_iter = 1e3;
    opti.solver('ipopt',options);

    % Solve problem
    sol = opti.solve();
    
    %% retrieve solution
    out.F    = sol.value(Frel);
    
    % copy
    newparms = parms;
    
    % parameters
    newparms.kappa      = sol.value(kappa);
    newparms.act_max    = sol.value(act_max);
    newparms.vmax       = sol.value(vmax);
    newparms.n          = sol.value(n);
    newparms.kpe        = sol.value(kpe);
    newparms.Lce0       = sol.value(Lce0);
    
    %% "simulate" output
    a = newparms.actfunc(Cas, newparms);
    
    Fce_sim = newparms.vF_func(vis * newparms.gamma, newparms) .* a;
    Fpe_sim = newparms.kpe * log(1+exp((Lts - newparms.Lce0)));

    Fse_sim = (Fce_sim + Fpe_sim) * newparms.Fscale;
    
    %% plot
    figure(2)
    plot(tis(idF), Fis(idF), '.', tis, out.F, '-', tis, Fse_sim, '--')
    
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

