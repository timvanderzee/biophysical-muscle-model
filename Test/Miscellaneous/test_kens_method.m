clear all; close all; clc

[username, githubfolder] = get_paths();

mcodes = [1 2 1; 1 2 1];
discretized_model = 1;
FLs = [0 1];

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

iFs = 7;
pCa = 6.1;
Ca = 10.^(-pCa+6);
AMP = 0;
ISI = .001;
parms_version = '_v3';

figure(1)
for iF = iFs
    nexttile
    
    
    % load parameters
    mcode = [1 2 1];
    [output_mainfolder, modelname, ~, ~] = get_folder_and_model(mcode);
    
    input_foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iF}];
    cd(input_foldername)
    load(['parms_',modelname, parms_version, '.mat'], 'newparms')
    parms = update_parms(newparms);
%     parms.FL_overlap = FLs(i);
    
    n0 = zeros(size(parms.xi));
    x0 = [n0'; parms.x0(4:end)'];
    xp0 = zeros(size(x0));
    X0 = x0;
    
    dTt = .0383/.4545; % test stretch (= constant)
    dTc = AMP / .4545; % conditioning stretch
    
    odeopt = odeset('maxstep', 1e-3);
    %         odeopt = [];
    %
    tiso = dTt*3+dTc*2+ISI + 5;
    dt = .001; % gives 10 points in SRS zone
    N = round(tiso / dt);
    
    [tis, Cas, Lis, vis, ts, Ts] = create_input(tiso, dTt, dTc, ISI, Ca, N);
    
    parms.ti = tis;
    parms.vts = vis;
    parms.Cas = mean(Cas);
    parms.Lts = Lis * parms.gamma;
    
    
    % run simulation
    
    aTs = [0; Ts];
    X0 = x0;
    vts = [0 .4545 -.4545 0 .4545 0 0];
    
    
    % interval needs to have finite duration
    nzi = find(diff(aTs) > 0);
    
    
    %% method of characteristics
    
    parms.xi = linspace(-15,15,1000);
    n0 = zeros(size(parms.xi));
    X0 = [n0'; parms.x0(4:end)'];

    tic
    
    Y = [];
    t1 = [];
    
    % method of characteristics
    parms.method_of_characteristics = 1;

    for p = 1:(length(nzi)-1)
        disp(p)        


        sol = ode45(@(t,y,yp) fiber_dynamics_explicit_no_tendon_full(t,y, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, odeopt);
        
        Y = [Y sol.y];
        t1 = [t1 sol.x];
        
        X0 = sol.y(:,end);
     
        Ns(p) = sol.stats.nfevals;

    end
    
     tot_time(1) = toc

     time_per_it(1) = tot_time(1) / sum(Ns);
    
    % compute force
    L = Y(end-3,:);
    F1 = nan(1, size(Y,2));
    n = Y(1:end-4,:);
    

    for i = 1:size(Y,2)
        xi = parms.xi + L(i);
        F1(i) = trapz(xi, xi .* n(:,i)') + trapz(xi, n(:,i)');
    end
        
    %% Ken's method
    % method of characteristics
    parms.method_of_characteristics = 0;
    
    dt = 1e-4;
    
    tlin = 0:dt:aTs(end);
    
    parms.xi = -3:.03:3;
    
    n0 = zeros(size(parms.xi));
    X0 = [n0'; parms.x0(4:end)'];
    
    y = X0;
    t2 = 0;
    
    N = length(tlin);
    F2 = zeros(1,N);
    
    vts = interp1(tis, vis, tlin);
    
    tic
    
    for id = 2:N
        
%         disp(tlin(id))
        
        parms.F = F2(id-1);
        parms.vts = vts(id);
        
        
        yp = fiber_dynamics_explicit_no_tendon_full(t2(id-1),y(:,id-1), parms);
        
        %%
%         profile on
        y(:,id) = y(:,id-1) + yp * dt;
        
        % for n, also need to shift
        dx = yp(end-3) * dt;
        xi = parms.xi - dx;
        n = y(1:length(parms.xi),id);
        
        nshift = interp1(parms.xi(:), n, xi(:));
        nshift(isnan(nshift)) = 0;
        nshift(nshift<0) = 0;
        
        y(1:length(parms.xi),id) = nshift;
        
        t2(id) = t2(id-1) + dt;
        
        F2(id) = trapz(parms.xi, parms.xi .* nshift') + trapz(parms.xi, nshift');
        
%         profile viewer
    end
    
    tot_time(2) = toc;
    
    time_per_it(2) = tot_time(2) / N;
    
    %% only do function evaluations
    tfval = [];
    parms.method_of_characteristics = 1;
    
    for id = 2:N
        
        tic
        yp = fiber_dynamics_explicit_no_tendon_full(t2(id-1),y(:,id-1), parms);
        tfval(id-1) = toc;
    end
    
    
    %%
    close all
    figure(1)
    plot(t1, F1, '-', t2, F2, '--')
    
end

% find unique values
[~, ui] = unique(tall);

% interpolate force
oFi = interp1(tall(ui), Fall(ui), tis) * parms.Fscale + parms.Fpe_func(parms.Lts, parms);



figure(1)
plot(tis, oFi, 'linewidth', 2); hold on
box off
xlabel('Time (s)')
ylabel('Force (-)')

%         xlim([4.5 max(tall)])



