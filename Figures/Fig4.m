clear all; close all; clc
usernames = {'timvd','u0167448'};

for i = 1:length(usernames)
    username = usernames{i};
    docfolder =  ['C:\Users\', usernames{i}, '\Documents'];
    
    if isfolder(docfolder)
        mainfolder = docfolder;
    end
end

if contains(mainfolder, 'timvd')
    githubfolder = mainfolder;
else
    githubfolder = [mainfolder, '\GitHub'];
end

addpath(genpath([githubfolder, '\muscle-thixotropy']))
addpath(genpath([mainfolder, '\casadi-3.7.1-windows64-matlab2018b']))
addpath(genpath([githubfolder, '\biophysical-muscle-model']))

figure(1)
color = get(gca,'colororder');
pcolors = flip(parula(7));
acolors = [color(2,:); pcolors(4:end-1,:);pcolors(4:end-1,:)];

fig = 6;

if fig == 4 || fig == 5
    mcodes = [2 1 1; 1 1 3; 1 1 1];
    vs = {'\', '\','\'};
    colors = acolors;

else
    mcodes = [1 1 1; 1 2 1];
    vs = {'\', '\full\'};
    colors = acolors(3:end,:);
end

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
if fig == 4
    ISIs = [.001 .001;
        .100 .100];
    
    AMPs = [0      0;
        .0383 .0383];
    
    pCas = [4.5 6.2;
        4.5 6.2];
    
    titles = {'Maximal activation', 'Submaximal activation'};
    
elseif fig == 5
    
    ISIs = [.001 .001;
        .001 .316];
    
    AMPs = [0      0;
        .0383 .0383];
    
    pCas = [6.2 6.2;
        6.2 6.2];
    
    
    titles = {'No recovery', 'Long recovery'};
    
elseif fig == 6
    
    ISIs = [.001 .001;
        .1 .316];
    
    AMPs = [0      0;
        .0383 .0383];
    
    pCas = [6.2 6.2;
        6.2 6.2];
    
    
    titles = {'Short recovery', 'Long recovery'};
end

for kk = 1:size(mcodes,1)
    disp(filenames{kk})
    output_folder = [opt_types{kk},'\normalized\with_PE_optimized\2_trials'];
    
    output_dir = [output_mainfolder{1}, '\', filenames{kk},vs{kk}, output_folder];
    cd(output_dir)
    
    load([filenames{kk},'_F', num2str(k),'_best.mat'],'parms','exitflag','fopt','C0','Cbounds','model','P0','P')
    C = p_to_c(P, Cbounds);
    parms = C_to_parms(C, parms, parms.optvars);
    parms = calc_dependent_parms(parms);
    
    parms.act = 1;
    parms.Noverlap = 1;
    parms.vF_func =  @(vcerel,parms)parms.e(1)*log((parms.e(2)*vcerel./parms.vmax+parms.e(3))+sqrt((parms.e(2)*vcerel./parms.vmax+parms.e(3)).^2+1))+parms.e(4);
    
    Parms{kk} = parms;
    
    %     if kk == 1
    %     cd([mainfolder, '\GitHub\biophysical-muscle-model\Parameters'])
    %     load('parms_v5.mat')
    %     Parms{kk} = sparms(k);
    %     end
end

%% evaluate
odeopt = odeset('maxstep', 1e-2);
gamma = 108.3333; % length scaling

ls = {'--','-'};

for j = 1:size(ISIs,1)
    %     if ishandle(j), close(j); end
    %     figure(j)
    
    figure(1)
    [texp, Lexp, Fexp, Tsrel] = get_data(data, ISIs(j,:), AMPs(j,:), pCas(j,:));
    plot_data(texp, Lexp, Fexp, brighten([.5 .5 .5], (2-j)/2), ls{j});
    
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
        
        %         subplot(3,3,i)
        %         plot(t(t<.15), Lis(t<.15)*100, 'linewidth',2)
        
        tids = sort([Tsrel(i,:) .15]);
        
        for kk = 1:size(mcodes,1)
            
            parms = Parms{kk};
            
            if contains(filenames{kk}, 'Hill')
                x0 = 0;
            else
                x0 = [parms.x0(2:end)'; 0];
            end
            
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
            
           
            figure(1)
            subplot(4,size(ISIs,2),[i+size(ISIs,2) i+size(ISIs,2)*3])
            plot(t(t<.15), oFi(t<.15)*100, 'linestyle', ls{j}, 'linewidth',2, 'color', brighten(colors(kk,:), (2-j)/2))
            
            % compute RMSD
            id = t < .15 & t > (-ISI - 2 * dTc - .1);
            oFii = interp1(t(id), oFi(id), texp(:,i));
            
            RMSDs = sqrt((oFii - Fexp(:,i)).^2) * 100;
            xlabel('Time (s)')
            
            figure(10)
            subplot(1,size(ISIs,2),i)
            plot(texp(:,i), RMSDs, 'linestyle', ls{j}, 'linewidth', 2, 'color', brighten(colors(kk,:), (2-j)/2));
            hold on
            box off
            
            RMSD = sqrt(sum((oFii - Fexp(:,i)).^2));
  
            for ii = 1:(length(tids)-2)
                plot([tids(ii) tids(ii)], [0 20], 'k:'); hold on
            end
            
            axis([-.3 .15 0 20])
            xlabel('Time (s)')
            
            if i == 1
                ylabel('RMSD (%F_0)')
            end
        end
        
        % add vertical lines
        figure(1)
        subplot(4,size(ISIs,2),i)
        title(titles{i})
        for ii = 1:(length(tids)-2)
            plot([tids(ii) tids(ii)], [0 20], 'k:'); hold on
        end
        
        subplot(4,size(ISIs,2),[i+size(ISIs,2) i+size(ISIs,2)*3])        
        for ii = 1:(length(tids)-2)
            plot([tids(ii) tids(ii)], [0 300], 'k:'); hold on
        end
        
        
    end
end

%%
figure(1)
subplot(4,size(ISIs,2),[1+size(ISIs,2) 1+size(ISIs,2)*3])       

if fig == 4 || fig == 5
legend('Data','Hill','XB','XB coop','location','best')
else
legend('Data','XB coop','XB coop + FD','location','best')
end
legend box off

%% optionally export to PNG
cd(['C:\Users\',username,'\OneDrive\9. Short-range stiffness\figures\MAT'])
       
figure(1)
exportgraphics(gcf,['Fig',num2str(fig),'.png'])









