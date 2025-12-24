clear all; close all; clc
savefig = 1;

[username, githubfolder] = get_paths();
th = [0 .07 .7 1.5];
tid = [3 1 7];

%% Data
Cid = tid(1);
ISIid = tid(2);
AMPid = tid(3);

% [username, githubfolder] = get_paths();
% cd([githubfolder, '\biophysical-muscle-model\Data'])
load('SRS_data_v3.mat', 'SRS_pre', 'F0', 'SRS_post')

% average some pCas
iFs = 1:11;
%  iFs = [1 2 3, 5, 6, 7, 8, 10, 11];
SRSrel = nan(length(th)-1,7,8,length(iFs));
F0s = nan(length(th)-1, 7,8,length(iFs));

for k = iFs
    for i = 1:length(th)-1
        id = F0(:,1,7,k) > th(i) & F0(:,1,7,k) <= th(i+1);
        
        %         SRSrel(i,:,:,k) = mean(SRS_post(id,:,:,k),1,'omitnan') ./ mean(SRS_pre(id,:,7,k),'all', 'omitnan');
        SRSrel(i,:,:,k) = mean(SRS_post(id,:,:,k),1,'omitnan') ./ mean(SRS_post(id,1,1,k),'all', 'omitnan');
        F0s(i,:,:,k) = mean(F0(id,:,:,k), 1,'omitnan');
    end
end

%%
% summary plot
eAMPs = [0 12 38 121 216 288 383 682]/10000;
eISIs = [1 10 100 316 1000 3160 10000]/1000;
aeAMPs = repmat(eAMPs, 7, 1);
aeISIs = repmat(eISIs(:), 1, 8);

% model
modelnames = {'Hill_alternative', 'Hill_regular', 'biophysical_full_alternative'};
titles = {'Hill model (no SE)', 'Hill model (with SE)', '4-state XB coop model'};
versions = {'parms', 'parms_v3', 'parms_v4d'};
% versions = {'parms_v3', 'parms_v4', 'parms_v6', 'parms_v4'};

% modelnames =  {'biophysical_full_alternative'};
acolors = lines(6);
acolors = [acolors(end,:); acolors];

id = [1, 2, 6];
color = acolors(id,:);

for ii = 1:length(modelnames)
    modelname = modelnames{ii};
    
    cd([githubfolder, '\biophysical-muscle-model\Model output\SRS\', versions{ii}]);
    load([modelname,'_SRS.mat'],'Stest', 'Scond', 'AMPs', 'iFs', 'pCas', 'ISIs', 'F0')
    
    amAMPs = repmat(AMPs, length(ISIs), 1);
    amISIs = repmat(ISIs(:), 1, length(AMPs));
    
    % average
    % iFsm = 9:11;
    iFsm = [2,3,5,6,7,8,11];
    iFsd = 1:11;
    
    SRSrel_m = nan(length(th)-1,length(ISIs),length(AMPs),iFs(end));
    F0s_m = nan(length(th)-1, length(ISIs),length(AMPs),iFs(end));
    
    miid = ismember(ISIs,  eISIs);
    maid = ismember(AMPs,  eAMPs);
    
    for k = iFsm
        for i = 1:length(th)-1
            id = F0(:,1,7,k) > th(i) & F0(:,1,7,k) <= th(i+1);
            
            %         SRSrel_m(i,:,:,k) = mean(Stest(id,:,:,k),1,'omitnan') ./ mean(Scond(id,1,7,k),'all', 'omitnan');
            SRSrel_m(i,:,:,k) = mean(Stest(id,:,:,k),1,'omitnan') ./ mean(Stest(id,1,1,k),'all', 'omitnan');
            F0s_m(i,:,:,k) = mean(F0(id,:,:,k), 1,'omitnan');
        end
    end
    
    % activation
    k = 2;
    
    fig = figure(1)
    subplot(3,3,ii)
    surf(amAMPs, amISIs, squeeze(mean(SRSrel_m(k,:,:,iFsm), 4, 'omitnan')),'edgecolor',[.8 .8 .8],'linestyle','-','Facealpha',.7,'linewidth',.5); hold on
    surf(aeAMPs,aeISIs, 0*ones(size(squeeze(aeISIs))),'facecolor', [1 1 1],'edgecolor', 'none'); hold on
    surf(aeAMPs(1:4,7:8),aeISIs(1:4,7:8), 0*ones(size(squeeze(aeISIs(1:4,7:8)))),'facecolor', [.9 .9 .9],'edgecolor', 'none'); hold on
    
    colors = color(ii,:) .* repmat(linspace(0,1,100)', 1, 3);
    set(gca, 'colormap', colors)
    
    plot3(aeAMPs(:,[5:6, 8]), aeISIs(:,[5:6, 8]), 0*ones(size(aeAMPs(:,[5:6, 8]))), 'color', [.8 .8 .8],'linewidth',2)
    plot3(aeAMPs([2 6],:)', aeISIs([2 6],:)', 0*ones(size(aeAMPs([2 6],:)))', 'color', [.8 .8 .8],'linewidth',2)
    
    plot3(aeAMPs(:,[2:4, 7]), aeISIs(:,[2:4, 7]), 0*ones(size(aeAMPs(:,[2:4, 7]))), 'color', [.6 .6 .6],'linewidth',2)
    plot3(aeAMPs([1 3:5 7],:)', aeISIs([1 3:5 7],:)', 0*ones(size(aeAMPs([1 3:5 7],:)))', 'color', [.6 .6 .6],'linewidth',2)
    
    plot3([aeAMPs(3,7) aeAMPs(3,7)], [aeISIs(3,7) aeISIs(3,7)], [0 squeeze(mean(SRSrel(k,3,7,iFsd), 4, 'omitnan'))], '-', 'linewidth',1.5, 'color', [.6 .6 .6])
    plot3([aeAMPs(1,1) aeAMPs(1,1)], [aeISIs(1,1) aeISIs(1,1)], [0 squeeze(mean(SRSrel(k,1,1,iFsd), 4, 'omitnan'))], '-', 'linewidth',1.5, 'color', [.6 .6 .6])
    
    set(gca,'YScale','log', 'Clim', [.6 1.1], 'xtick', 0:.02:.1, 'ytick', [1e-3 1e-1 1e1]);
    zlim([0 1.1])
    axis([0 .07 8e-4 1.5e1 0 1.1])
    
    plot3(aeAMPs, aeISIs, squeeze(mean(SRSrel(k,:,:,iFsd), 4, 'omitnan')),'o', 'color', [.5 .5 .5], 'markerfacecolor', [.5 .5 .5],'markersize',3); hold on
    
    plot3(amAMPs(1,:), amISIs(1,:), squeeze(mean(SRSrel_m(k,1,:,iFsm), 4, 'omitnan')),'-','linewidth',1.5, 'color', color(ii,:))
    plot3(amAMPs(:,1), amISIs(:,1), squeeze(mean(SRSrel_m(k,:,1,iFsm), 4, 'omitnan')),'-','linewidth',1.5, 'color', color(ii,:))
    plot3(amAMPs(end,:), amISIs(end,:), squeeze(mean(SRSrel_m(k,end,:,iFsm), 4, 'omitnan')),'-','linewidth',1.5, 'color', color(ii,:))
    plot3(amAMPs(:,end), amISIs(:,end), squeeze(mean(SRSrel_m(k,:,end,iFsm), 4, 'omitnan')),'-','linewidth',1.5, 'color', color(ii,:))
    
    set(gca,'YScale','log')
    
    %% compute R2
    SSE = sum((mean(SRSrel_m(k,miid,maid,iFsm), 4, 'omitnan') - mean(SRSrel(k,:,:,iFsd), 4, 'omitnan')).^2,'all', 'omitnan')
    SST = sum((mean(SRSrel(k,:,:,iFsd), 'all','omitnan') - mean(SRSrel(k,:,:,iFsd), 4, 'omitnan')).^2,'all', 'omitnan')
    
    R2 = 1 - SSE./SST
    
    set(gca,'fontsize',6)
    title(titles{ii}, 'fontsize', 8)
    if ii == 1
        zlabel('Relative short-range stiffness','Fontsize',8)
    else
        set(gca, 'zticklabel', {}, 'ZColor', 'none')
    end
    
    view(230-180,20)
    
    text(.07, 1e0, 1.4, ['R^2 = ',num2str(round(R2, 2), 3)], 'fontsize', 6, 'horizontalalignment', 'center')
    
    
end

%%
figure(1)

subplot(331)
h = axes(fig,'visible','off');
set(h, 'Clim', [.6 1.1])
% g = colorbar(h, 'location', 'WestOutside')

% set(g,'units','centimeters','position',[6 2.5+11 0.2 1.7], 'fontsize', 6)

%% plot force traces
% close all

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
% cd([output_mainfolder{2},'\data'])

% if ishandle(2), close(2); end
% figure(2)
iFsd = 6;


mcodes = [2 2 1; 2 1 1; 1 1 1];
versions = {'parms', 'parms_v3', 'parms_v4d', 'parms_v4d'};

pCas = [9 6.1 4.5];

dy = 0;

ls = {'--',':','-'};

ax(1) = subplot(3,3,4);
ax(2) = subplot(3,3,5);
ax(3) = subplot(3,3,6);
ax(4) = subplot(3,3,[7 9]);

aISIs = [.001 .001 .001 1];
AMPs = .0383;

for k = 1:4
    ISIs = aISIs(k);
    set(ax(k), 'fontsize', 6);
    
    if k == 1 || k == 4
        pCas =  [9 6.1 4.5];
        mcodes = [2 2 1; 2 1 1; 1 1 1];
        versions = {'parms', 'parms_v3', 'parms_v4d', 'parms_v4d'};
        ls = {'--',':','-'};
        id = [1, 2, 6];
        color = acolors(id,:);

    else
        pCas = 6.1;
        mcodes = [2 1 1; 1 1 1];
        versions = {'parms_v3', 'parms_v4d', 'parms_v4d'};
        ls = {':','-'};
        id = [2, 6];
        color = acolors(id,:);
    end
    
    for i = 1:length(pCas)
        pCa = pCas(i);
        
        cd(['C:\Users\',username,'\OneDrive - KU Leuven\9. Short-range stiffness\matlab\data'])
        %         load([fibers{iF},'_cor_new.mat'],'data')
        

        for iF = iFsd
            load([fibers{iF},'_cor_new.mat'],'data');
            [texp, Lexp, Fexp, Tsrel] = get_data(data, ISIs, AMPs, pCa);
            
            plot(ax(k), texp - .003, Fexp + dy * i,'k-', 'linewidth', 2);
            hold(ax(k), 'on')
        end
        
        
        for j = flip(1:size(mcodes,1))
            mcode = mcodes(j,:);
            
            [output_mainfolder, modelfilename, ~, ~] = get_folder_and_model(mcode);
            dTt = .0383/.4545; % test stretch (= constant)
            dTc = AMPs / .4545; % conditioning stretch
            ISI = ISIs;
            AMP = AMPs;
            
            %
            %             clc
            tiso = dTt*3+dTc*2+ISI + 2;
            
            
            filename = [output_mainfolder{2}, '\', versions{j},'\', modelfilename,'\',fibers{iF}, '\pCa=',num2str(pCa*10),'\', fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'];
            disp(filename)
            
            
            load(filename, 'tis','Cas','vis','Lis','oFi','parms', 'ts')
            
            t = tis + 3*dTt - tiso;
            plot(ax(k), t, oFi + dy * i, 'linewidth', 2, 'color', color(j,:), 'linestyle', ls{j})
            
            xlim(ax(k), [-.2-ISIs .15])
        end
        box(ax(k), 'off')
        %         ylim([0 3])
        
    end
end

%% make nice

figure(1)

set(ax(3), 'ycolor', [0 0 0])



for k = 1:4
    
    ax(k).YLabel.String = 'Force (-)';
    ax(k).YLabel.Color = [0 0 0];
    ax(k).YLabel.FontSize = 8;
    
    ax(k).XLabel.String = 'Time (s)';
    ax(k).XLabel.Color = [0 0 0];
    ax(k).XLabel.FontSize = 8;
    
end

title(ax(2), 'Before shortening', 'fontsize', 8)
title(ax(1), 'Trial with short recovery', 'fontsize', 8)
title(ax(3), 'After shortening', 'fontsize', 8)
title(ax(4), 'Trial with long recovery', 'fontsize', 8)

%% zoom in
P1 = [-.175 -.125 .5 1];
P2 = [.05 .1 .5 1];

figure(1)
axis(ax(1), [-.2 .12 0 2])
axis(ax(3), P2)
axis(ax(2), P1)

h1 = rectangle(ax(1), 'position', [P1(1) P1(3) P1(2)-P1(1) P1(4)-P1(3)]);
h3 = rectangle(ax(1), 'position', [P2(1) P2(3) P2(2)-P2(1) P2(4)-P2(3)]);

h2 = rectangle(ax(2), 'position', [P1(1) P1(3) P1(2)-P1(1) P1(4)-P1(3)]);
h4 = rectangle(ax(3), 'position', [P2(1) P2(3) P2(2)-P2(1) P2(4)-P2(3)]);

yline(ax(1), .87, 'k-.')
yline(ax(2), .87, 'k-.')
yline(ax(3), .87, 'k-.')
yline(ax(4), .87, 'k-.')

%% size
figure(1)
set(gcf,'units','centimeters','position',[10 2 13 20])

%%
if savefig
    cd(['C:\Users\',username,'\OneDrive\9. Short-range stiffness\figures\MAT'])
    figure(1)
    exportgraphics(gcf,['Fig9.png'])
end

