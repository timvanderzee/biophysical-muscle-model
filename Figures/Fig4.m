clear all; close all; clc
[username, githubfolder] = get_paths();

figure(1)
color = get(gca,'colororder');
pcolors = flip(parula(7));
acolors = [color(2,:); pcolors(4:end-1,:);pcolors(4:end-1,:)];

%% chose figure number: specify conditions
fig = 6;
iF = 6;

% chosen ISIs, AMPs and pCas
if fig == 4
    ISIs = [.001 .001; % dashed
            .100 .100]; % solid
    
    AMPs = [0      0; % dashed
        .0383 .0383]; % solid
    
    pCas = [4.5 6.2;
            4.5 6.2];
    
    titles = {'Maximal activation', 'Submaximal activation'};
    
elseif fig == 5
    
    ISIs = [.316 .316;
            .001 .001];
    
    AMPs = [.0383 .0383;
            .0383 .0383];
    
    pCas = [4.5 6.2;
            4.5 6.2];
    
    titles = {'Maximal activation', 'Submaximal activation'};
    
elseif fig == 6
    
    ISIs = [.001 .316;
            .100 .001]; % solid
    
    AMPs = [0   .0383;
            .0383       .0383];
    
    pCas = [6.2 6.2;
            6.2 6.2];
    
    titles = {'Effect of amplitude', 'Effect of recovery'};
end

ISIs = flip(ISIs,1);
AMPs = flip(AMPs,1);
pCas = flip(pCas,1);

%% choose fiber: load data and parameters
if fig == 4 || fig == 5
    mcodes = [2 1 1; 1 1 3; 1 1 1];
    colors = acolors;

else
    mcodes = [1 1 1; 1 2 1];
    colors = acolors(3:end,:);
end

for kk = 1:size(mcodes,1)
    mcode = mcodes(kk,:);
    
    [output_mainfolder, filenames{kk}, opt_types{kk}, ~] = get_folder_and_model(mcode);

    cd([githubfolder, '\biophysical-muscle-model\Parameters'])
    load(['parms_',filenames{kk},'.mat'], 'pparms')
    Parms{kk} = pparms(iF);
end

% load data 
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
cd([output_mainfolder{2},'\data'])
load([fibers{iF},'_cor_new.mat'],'data');

%% evaluate
odeopt = odeset('maxstep', 1e-2);
gamma = 108.3333; % length scaling

ls = {'-',':'};

for j = 1:size(ISIs,1)
    %     if ishandle(j), close(j); end
    %     figure(j)
    
    figure(1)
    [texp, Lexp, Fexp, Tsrel] = get_data(data, ISIs(j,:), AMPs(j,:), pCas(j,:));
    plot_data(texp, Lexp, Fexp, brighten([.5 .5 .5], (j-1)/2), ls{j});
    
    for i = 1:size(ISIs,2)
        
        % evaluate fit
        Ca = 10.^(-pCas(j,i)+6);
        
        dTt = .0383/.4545; % test stretch (= constant)
        dTc = AMPs(j,i) / .4545; % conditioning stretch
        ISI = ISIs(j,i);
        AMP = AMPs(j,i);
%         tiso = 3;
        
%         [tis, Cas, Lis, vis, ts] = create_input(tiso, dTt, dTc, ISI, Ca, 2000);
        
        % center around 2nd stretch
%         t = tis + 3*dTt - tiso;
        
        %         subplot(3,3,i)
        %         plot(t(t<.15), Lis(t<.15)*100, 'linewidth',2)
        
        tids = sort([Tsrel(i,:) .15]);
        
        for kk = 1:size(mcodes,1)
            


%             ISI = ISIs(jj);

            tiso = dTt*3+dTc*2+ISI + 2;

            cd([output_mainfolder{2}])

            [output_mainfolder, filename, ~, ~] = get_folder_and_model(mcodes(kk,:));
            cd([filename,'\',fibers{iF}, '\pCa=',num2str(pCas(j,i)*10)])

            disp([fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'])
            load([fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'], ...
                'tis','Cas','vis','Lis','oFi','parms', 'ts')
            
            t = tis + 3*dTt - tiso;
%             parms = Parms{kk};
%             
%             if contains(filenames{kk}, 'Hill')
%                 x0 = 0;
%             else
%                 x0 = [parms.x0(2:end)'; 0];
%             end
%             
%             xp0 = zeros(size(x0));
%            
%             
%             parms.ti = tis;
%             parms.vts = vis;
%             parms.Cas = Cas;
%             parms.Lts = Lis;
%             
%             % run simulation
%             tic
%             if contains(filenames{kk}, 'Hill')
%                 sol = ode15i(@(t,y,yp) hill_type_implicit(t,y,yp, parms), [0 max(tis)], x0, xp0, odeopt);
%             else
%                 sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 max(tis)], x0, xp0, odeopt);
%             end
%             toc
%             
%             % update x0
%             x0 = sol.y(:,end);
%             
%             Liss = Lis * gamma;
%             ot = sol.x;
%             
%             if contains(filenames{kk}, 'Hill')
%                 lce = sol.y(1,:);
%                 
%                 % elastic elements
%                 Lse = Liss - interp1(sol.x, lce, tis);
%                 oFi = parms.Fse_func(Lse, parms) * parms.Fscale + parms.Fpe_func(Liss, parms);
%                 
%             else
%                 F = sol.y(1,:) + sol.y(2,:);
%                 
%                 oF = F * parms.Fscale;
%                 
%                 oFi = interp1(ot, oF, tis) + parms.Fpe_func(Liss, parms);
%             end
            
           
            figure(1)
            subplot(4,size(ISIs,2),[i+size(ISIs,2) i+size(ISIs,2)*3])
            plot(t(t<.15), oFi(t<.15)*100, 'linestyle', ls{j}, 'linewidth',2, 'color', brighten(colors(kk,:), (j-1)/2))
            
            % compute RMSD
            id = t < .15 & t > (-ISI - 2 * dTc - .1);
            oFii = interp1(t(id), oFi(id), texp(:,i));
            
            RMSDs = sqrt((oFii - Fexp(:,i)).^2) * 100;
            xlabel('Time (s)')
            
            figure(10)
            subplot(1,size(ISIs,2),i)
            plot(texp(:,i), RMSDs, 'linestyle', ls{j}, 'linewidth', 2, 'color', brighten(colors(kk,:), (j-1)/2));
            hold on
            box off
            
            RMSD = sqrt(sum((oFii - Fexp(:,i)).^2));
  
            for ii = 1:(length(tids)-2)
                plot([tids(ii) tids(ii)], [0 20], ':', 'color', [.5 .5 .5]); hold on
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
            plot([tids(ii) tids(ii)], [0 20], ':', 'color', [.5 .5 .5]); hold on
        end
        
        subplot(4,size(ISIs,2),[i+size(ISIs,2) i+size(ISIs,2)*3])        
        for ii = 1:(length(tids)-2)
            plot([tids(ii) tids(ii)], [0 300], ':', 'color', [.5 .5 .5]); hold on
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
% exportgraphics(gcf,['Fig',num2str(fig),'.png'])









