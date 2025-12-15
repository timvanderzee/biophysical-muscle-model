clear all; close all; clc
[username, githubfolder] = get_paths();
savefig = 0;

figure(1)
color = get(gca,'colororder');
pcolors = flip(parula(7));
acolors = lines(4);
% acolors = [color(2,:); pcolors(4:end-1,:);pcolors(4:end-1,:)];

discretized_model = 0;

%% chose figure number: specify conditions
fig = 6;
iF = 7;

pCas = [4.5 6.2 6.4;
        4.5 6.2 6.4];

% chosen ISIs, AMPs and pCas
if fig == 4
    ISIs = [.100 .100 .100;
            .001 .001 .001]; 
    
    AMPs = [.0383 .0383 .0383;
            0      0 0];% solid
    
    
    titles = {'Maximal activation', 'Submaximal activation'};
    
    tmax =  0.3195;
    
elseif fig == 5 || fig == 6
    
    ISIs = [ .001 .001 .001;
             .316 .316 .316];
    
    AMPs = [.0383 .0383 .0383;
            .0383 .0383 .0383];
        
        tmax = 0.6345;
    
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
    mcodes = [1 1 1; 1 1 1];
    colors = acolors(3:4,:);
    
    
    if discretized_model
        versions = {'_v2d', '_v2d'};
    else
%         versions = {'_v4', '_v4'};
        versions = {'_v4', '_v6'};
    end
end

% load data
[output_mainfolder, filename, ~, ~] = get_folder_and_model(mcodes(1,:));
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
cd([output_mainfolder{2},'\data'])
load([fibers{iF},'_cor_new.mat'],'data');

%% evaluate
odeopt = odeset('maxstep', 1e-2);
ls = {':','-'};

for j = 1:size(ISIs,1)
    %     if ishandle(j), close(j); end
    %     figure(j)
    
    figure(1)
    [texp, Lexp, Fexp, Tsrel] = get_data(data, ISIs(j,:), AMPs(j,:), pCas(j,:));
%     plot_data(texp, Lexp, Fexp, brighten([.5 .5 .5], (j-1)/2), ls{j});
    
    texp = texp - Tsrel(1,2);
%     tmax = .15 - Tsrel(1,2)

    for i = 1:size(ISIs,2)
    
        id = texp(:,i) < tmax;
        
        subplot(5,1,1)
        plot(texp(id,i), Lexp(id,i)*100, 'color', brighten([.5 .5 .5], (2-j)/2), 'linewidth', 2, 'linestyle',ls{j}); hold on
        axis([-.05 tmax -.5 5])
        box off
        xticks(-.05:.05:tmax);

        if j == 1
            set(gca, 'Fontsize', 6)
            ylabel('Length (%L_0)', 'Fontsize', 8)
        end
        
        subplot(5,1,[2 5])
        plot(texp(id,i), Fexp(id,i)*100, 'color', brighten([.5 .5 .5], (2-j)/2), 'linewidth', 2, 'linestyle',ls{j}); hold on
        axis([-.05 tmax 0 200])
        xticks(-.05:.05:tmax);
        
        box off

        if j == 1
            set(gca, 'Fontsize', 6)
            ylabel('Force (%F_0)', 'Fontsize', 8)
        end
        
        % evaluate fit
        Ca = 10.^(-pCas(j,i)+6);
        
        dTt = .0383/.4545; % test stretch (= constant)
        dTc = AMPs(j,i) / .4545; % conditioning stretch
        ISI = ISIs(j,i);
        AMP = AMPs(j,i);
        
        tids = sort([Tsrel(i,:) .15]) - Tsrel(1,2);
        
        for kk = 1:size(mcodes,1)
            mcode = mcodes(kk,:);
            [output_mainfolder, modelfilename, ~, ~] = get_folder_and_model(mcodes(kk,:));
            tiso = dTt*3+dTc*2+ISI + 2;
            
            filename = [output_mainfolder{2}, '\parms', versions{kk},'\', modelfilename,'\',fibers{iF}, '\pCa=',num2str(pCas(j,i)*10),'\', fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'];
            disp(filename)
            
            if ~exist(filename, 'file')
                disp('using approximated')
                filename = [output_mainfolder{2}, '\parms_v4\', modelfilename,'\',fibers{iF}, '\pCa=',num2str(pCas(j,i)*10),'\', fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'];
                
            end
            
            load(filename, 'tis','Cas','vis','Lis','oFi','parms', 'ts')
            
            t = tis + 3*dTt - tiso;
            
            t = t - Tsrel(1,2);
            tmax2 = .15 - Tsrel(1,2);
            
            figure(1)
            subplot(5,1,[2 5])
            plot(t(t<tmax2), oFi(t<tmax2)*100, 'linestyle', ls{j}, 'linewidth',2, 'color', brighten(colors(kk,:), (2-j)/2)); hold on
            
            % compute RMSD
            id = t < .15 & t > (-ISI - 2 * dTc - .1);
            oFii = interp1(t(id), oFi(id), texp(:,i));
            
            xlabel('Time (s)', 'Fontsize', 8)

        end
        
        % add vertical lines
        if ISI > .01
            phases = {'Isometric','Pre-stretch','Shortening','Recovery','Test stretch','Isometric'};
        else
            phases = {'Isometric','Pre-stretch','Shortening','','Test stretch','Isometric'};
        end
        
        figure(1)
        subplot(5,1,1)
%         title(titles{1}, 'Fontsize', 8)
        for ii = 1:(length(tids)-2)
            plot([tids(ii) tids(ii)], [0 20], ':', 'color', [.5 .5 .5]); hold on
            
            if j == 2
                text(mean([tids(ii) tids(ii+1)]), 4.5, phases{ii}, 'Fontsize',8,'HorizontalAlignment','center')
            end
        end
        
        subplot(5,1,[2 5])
        for ii = 1:(length(tids)-2)
            xline([tids(ii) tids(ii)], ':', 'color', [.5 .5 .5]); hold on
        end
        
        
    end
end

%%
figure(1)
subplot(5,1,[2 5])

if fig == 4
    legend('Data','Hill','2-state','3-state coop','location','north', 'Fontsize', 8)
elseif fig == 5
    legend('Data','Hill','2-state','3-state coop','location','best', 'Fontsize', 8)
else
    legend('Data','3-state coop','4-state coop','location','best', 'Fontsize', 8)
end
% legend box off

% subplot(4,1,1)
% text(-.4, 6, 'A', 'fontsize', 12, 'fontweight', 'bold')

% subplot(4,2,2)
% text(-.4, 6, 'B', 'fontsize', 12, 'fontweight', 'bold')

%% size
figure(1)

trange = [-.05 tmax];


set(gcf,'units','centimeters','position',[10 5 diff(trange)*25 15])
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







