clear all; close all; clc
[username, githubfolder] = get_paths();
savefig = 1;

figure(1)
color = get(gca,'colororder');
pcolors = flip(parula(7));
acolors = [color(2,:); pcolors(4:end-1,:);pcolors(4:end-1,:)];

discretized_model = 0;

%% chose figure number: specify conditions
fig = 6;
iF = 6;

pCas = [4.5 6.2;
    4.5 6.2];

% chosen ISIs, AMPs and pCas
if fig == 4
    ISIs = [.001 .001; % dashed
        .100 .100]; % solid
    
    AMPs = [0      0; % dashed
        .0383 .0383]; % solid
    
    
    titles = {'Maximal activation', 'Submaximal activation'};
    
elseif fig == 5 || fig == 6
    
    ISIs = [.316 .316;
        .001 .001];
    
    AMPs = [.0383 .0383;
        .0383 .0383];
    
    titles = {'Maximal activation', 'Submaximal activation'};
    
end

ISIs = flip(ISIs,1);
AMPs = flip(AMPs,1);
pCas = flip(pCas,1);

%% choose fiber: load data and parameters
if fig == 4 || fig == 5
    mcodes = [2 1 1; 1 1 3; 1 1 1];
    colors = acolors;
    
    %     versions = {'_v3', '_v3', '_v4'};
    
    if discretized_model
        versions = {'_v3', '_v2d', '_v2d'};
    else
        versions = {'_v3', '_v3', '_v4'};
    end
    
else
    mcodes = [1 1 1; 1 2 1];
    colors = acolors(3:end,:);
    
    
    if discretized_model
        versions = {'_v2d', '_v2d'};
    else
        versions = {'_v4', '_v4'};
    end
end

% load data
[output_mainfolder, filename, ~, ~] = get_folder_and_model(mcodes(1,:));
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
cd([output_mainfolder{2},'\data'])
load([fibers{iF},'_cor_new.mat'],'data');

%% evaluate
odeopt = odeset('maxstep', 1e-2);
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
        
        tids = sort([Tsrel(i,:) .15]);
        
        for kk = 1:size(mcodes,1)
            mcode = mcodes(kk,:);
            [output_mainfolder, modelfilename, ~, ~] = get_folder_and_model(mcodes(kk,:));
            tiso = dTt*3+dTc*2+ISI + 2;
            
            %             cd([output_mainfolder{2}])
            %             cd(['parms', versions{kk}])
            
            
            %             [output_mainfolder, filename, ~, ~] = get_folder_and_model(mcodes(kk,:));
            %             cd([output_mainfolder{2}][filename,'\',fibers{iF}, '\pCa=',num2str(pCas(j,i)*10)])
            filename = [output_mainfolder{2}, '\parms', versions{kk},'\', modelfilename,'\',fibers{iF}, '\pCa=',num2str(pCas(j,i)*10),'\', fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'];
            disp(filename)
            
            if ~exist(filename, 'file')
                disp('using approximated')
                filename = [output_mainfolder{2}, '\parms_v4\', modelfilename,'\',fibers{iF}, '\pCa=',num2str(pCas(j,i)*10),'\', fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'];
                
            end
            
            load(filename, 'tis','Cas','vis','Lis','oFi','parms', 'ts')
            
            t = tis + 3*dTt - tiso;
            
            figure(1)
            subplot(4,size(ISIs,2),[i+size(ISIs,2) i+size(ISIs,2)*3])
            plot(t(t<.15), oFi(t<.15)*100, 'linestyle', ls{j}, 'linewidth',2, 'color', brighten(colors(kk,:), (j-1)/2))
            
            % compute RMSD
            id = t < .15 & t > (-ISI - 2 * dTc - .1);
            oFii = interp1(t(id), oFi(id), texp(:,i));
            
            %             set(gca, 'Fontsize', 6)
            
            %             RMSDs = sqrt((oFii - Fexp(:,i)).^2) * 100;
            xlabel('Time (s)', 'Fontsize', 8)
            
            %             figure(10)
            %             subplot(1,size(ISIs,2),i)
            %             plot(texp(:,i), RMSDs, 'linestyle', ls{j}, 'linewidth', 2, 'color', brighten(colors(kk,:), (j-1)/2));
            %             hold on
            %             box off
            %
            %             RMSD = sqrt(sum((oFii - Fexp(:,i)).^2));
            %
            %             for ii = 1:(length(tids)-2)
            %                 plot([tids(ii) tids(ii)], [0 20], ':', 'color', [.5 .5 .5]); hold on
            %             end
            %
            %             axis([-.3 .15 0 20])
            %             xlabel('Time (s)')
            %
            %             if i == 1
            %                 ylabel('RMSD (%F_0)')
            %             end
        end
        
        % add vertical lines
        if ISI > .01
            phases = {'Isometric','Pre-stretch','Shortening','Recovery','Test stretch','Isometric'};
        else
            phases = {'Isometric','Pre-stretch','Shortening','','Test stretch','Isometric'};
        end
        
        figure(1)
        subplot(4,size(ISIs,2),i)
        title(titles{i}, 'Fontsize', 8)
        for ii = 1:(length(tids)-2)
            plot([tids(ii) tids(ii)], [0 20], ':', 'color', [.5 .5 .5]); hold on
            
            if j == 1
                text(mean([tids(ii) tids(ii+1)]), 4.5, phases{ii}, 'Fontsize',4,'HorizontalAlignment','center')
            end
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

if fig == 4
    legend('Data','Hill','XB','XB coop','location','Southwest', 'Fontsize', 8)
elseif fig == 5
    legend('Data','Hill','XB','XB coop','location','South', 'Fontsize', 8)
else
    legend('Data','XB coop','XB coop + FD','location','South', 'Fontsize', 8)
end
legend box off

subplot(4,size(ISIs,2),1)
text(-.4, 6, 'A', 'fontsize', 12, 'fontweight', 'bold')

subplot(4,size(ISIs,2),2)
text(-.4, 6, 'B', 'fontsize', 12, 'fontweight', 'bold')

%% size
figure(1)

set(gcf,'units','centimeters','position',[10 10 19 10])
% set(gcf,'units','centimeters','position',[5 2.5 0.2 1.7], 'fontsize', 6)

% pause(0.5)

%% optionally export to PNG
cd(['C:\Users\',username,'\OneDrive\9. Short-range stiffness\figures\MAT'])

if discretized_model
    figname = ['Fig',num2str(fig),'.png'];
else
    figname = ['FigS',num2str(fig-3),'.png'];
end


if savefig
    figure(1)
    exportgraphics(gcf,figname)
end







