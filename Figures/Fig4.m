clear all; close all; clc

% mainfolder = 'C:\Users\u0167448\Documents';
mainfolder = 'C:\Users\timvd\Documents';
username = 'timvd';
% username = 'u00167448';

if strcmp(username, 'timvd')
    addpath(genpath([mainfolder, '\muscle-thixotropy']))
    addpath(genpath([mainfolder, '\casadi-3.7.1-windows64-matlab2018b']))
    addpath(genpath([mainfolder, '\biophysical-muscle-model']))

else
    
    addpath(genpath([mainfolder, '\GitHub\muscle-thixotropy']))
    addpath(genpath([mainfolder, '\casadi-3.7.1-windows64-matlab2018b']))
    addpath(genpath([mainfolder, '\GitHub\biophysical-muscle-model']))
end

mcodes = [1 1 1; 1 1 3; 2 1 1];
      
for kk = 1:size(mcodes,1)
    mcode = mcodes(kk,:);
    
    [output_mainfolder, filenames{kk}, opt_types{kk}, ~] = get_folder_and_model(mcode);
end
%% load data and parameters
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

cd([output_mainfolder{2},'\data'])

k = 6;
load([fibers{k},'_cor_new.mat'],'data');

% chosen ISIs, AMPs and pCas
ISIs = [.001 .100 .316; 
        .001 .001 .001;
        .001 .001 .001];

AMPs = [.0383 .0383 .0383; 
        .0038 .0121 .0383; 
        .0383 .0383 .0383];

pCas = [6.1 6.1 6.1; 
        6.1 6.1 6.1; 
        6.2 6.1 4.5];

vs = {'\', '\'};

for kk = 1:size(mcodes,1)
    disp(filenames{kk})
    output_folder = [opt_types{kk},'\normalized\with_PE_optimized\2_trials'];

    output_dir = [output_mainfolder{1}, '\', filenames{kk},vs{1}, output_folder];
    cd(output_dir)

    load([filenames{kk},'_F', num2str(k),'_best.mat'],'parms','exitflag','fopt','C0','Cbounds','model','P0','P')
    C = p_to_c(P, Cbounds);
    parms = C_to_parms(C, parms, parms.optvars);
    parms = calc_dependent_parms(parms);  

    parms.act = 1;
    parms.Noverlap = 1;
    parms.vF_func =  @(vcerel,parms)parms.e(1)*log((parms.e(2)*vcerel./parms.vmax+parms.e(3))+sqrt((parms.e(2)*vcerel./parms.vmax+parms.e(3)).^2+1))+parms.e(4);

    Parms{kk} = parms;
end

%% evaluate
odeopt = odeset('maxstep', 1e-2);
gamma = 108.3333; % length scaling

for j = 1:size(ISIs,1)
    if ishandle(j), close(j); end
    figure(j)
    plot_data(data, ISIs(j,:), AMPs(j,:), pCas(j,:))

    for i = 1:size(ISIs,2)
        
        % evaluate fit
        Ca = 10.^(-pCas(j,i)+6);

        dTt = .0383/.4545; % test stretch (= constant)
        dTc = AMPs(j,i) / .4545; % conditioning stretch
        ISI = ISIs(j,i);
        tiso = 3;

        [tis, Cas, Lis, vis, ts] = create_input(tiso, dTt, dTc, ISI, Ca, 2000);
        
        % center around 2nd stretch
        t = tis + 3*dTt - tiso;
        
        subplot(2,3,i)
        plot(t(t<.15), Lis(t<.15)*100, 'linewidth',2)
        
        for kk = 1:size(mcodes,1)
        
            parms = Parms{kk};

            x0 = parms.x0(2:end)';
            xp0 = zeros(size(x0));


            parms.ti = tis;
            parms.vts = vis;
            parms.Cas = Cas;
            parms.Lts = Lis;

            % run simulation
            tic
            if contains(filenames{kk}, 'Hill')
                sol = ode15i(@(t,y,yp) hill_type_implicit(t,y,yp, parms), [0 max(tis)], x0, xp0, odeopt);
            else
                sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 max(tis)], x0, xp0, odeopt);
            end
            toc

            % update x0
            x0 = sol.y(:,end);

            Liss = Lis * gamma;
            ot = sol.x;

            if contains(filenames{kk}, 'Hill')
                lce = sol.y(1,:);

                % elastic elements
                Lse = Liss - interp1(sol.x, lce, tis);    
                oFi = parms.Fse_func(Lse, parms) * parms.Fscale + parms.Fpe_func(Liss, parms);

            else
                F = sol.y(1,:) + sol.y(2,:);

                oF = F * parms.Fscale;

                oFi = interp1(ot, oF, tis) + parms.Fpe_func(Liss, parms);
            end

            subplot(2,3,i+3)
            plot(t(t<.15), oFi(t<.15)*100, 'linewidth',2)
        end
        
    end
end

%%
function[] = plot_data(data, sISIs, sAMPs, spCas)

% all ISIs and AMPs
AMPs = [0 12 38 121 216 288 383 682]/10000;
ISIs = [1 10 100 316 1000 3160 10000]/1000;
pCas = [4.5000    6.1000    6.2000    6.3000    6.4000    6.6000    9.0000];

% velocity 
v = .4545; % L0/s
tiso = 3;

% XData = nan(2e4, length(sISIs));
Fexp = nan(2e4, length(sISIs));
texp = nan(2e4, length(sISIs));
Lexp = nan(2e4, length(sISIs));
Tsrel = nan(length(sISIs), 7);

for kk = 1:length(sISIs)
    ISI = sISIs(kk);
    AMP = sAMPs(kk);
    pCa = spCas(kk);
    
    n = find(ISIs == round(ISI,4));
    m = find(AMPs == round(AMP,4));
    i = find(pCas == round(pCa,1));

    if isempty(data.Fexp(:,i,n,m))
        keyboard
    end
    
     % selected data
    Fexp(:,kk) = data.Fexp(:,i,n,m);
    texp(:,kk) = data.texp(:,i,n,m);
    Lexp(:,kk) = data.Lexp(:,i,n,m);
    
    T = [AMP AMP .0383]./v; % time of stretch (s)
    Ts = [tiso T(1) T(2) ISI T(3) tiso];
    
    Tsrel(kk,:) = [-10 cumsum(Ts(1,:)) - sum(Ts(1,1:4))];

end


%% plot
for kk = 1:length(sISIs)
    tids = sort([Tsrel(kk,:) .15]);
    
    subplot(2,3,kk)
    plot(texp(:,kk), Lexp(:,kk)*100, 'color', [.5 .5 .5], 'linewidth', 2); hold on
    
    for ii = 1:(length(tids)-2)
        plot([tids(ii) tids(ii)], [0 200], 'k:'); hold on
    end
    
    axis([-.35 .25 0 5])
    box off
    
    if kk == 1
        ylabel('Length (%L_0)')
    end
    
    
    subplot(2,3,kk+3)
    plot(texp(:,kk), Fexp(:,kk)*100, 'color', [.5 .5 .5], 'linewidth', 2); hold on
    
    for ii = 1:(length(tids)-2)
        plot([tids(ii) tids(ii)], [0 200], 'k:'); hold on
    end
    
    axis([-.35 .25 0 max(Fexp(:,kk)*100)*1.2])
    box off
    
    xlabel('Time (s)')
    
    if kk == 1
        ylabel('Force (%F_0)')
    end
        
end
end

