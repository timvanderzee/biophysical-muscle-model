clear all; close all; clc
savefig = 1; 

[username, githubfolder] = get_paths();
mcodes = [2 1 1; 1 1 3; 1 1 1; 1 2 1];

% iFs = [1 2 3, 5, 6, 7, 8, 10, 11];

iFs = [2,3,5,6,7,8,11];

th = [0 .07 .25 .7 1.5];

RMSDs = nan(length(th)-1,7,8,11, 7, 4);
F0s = nan(length(th)-1, 7,8,11, 4);
RMSDc = nan(length(th)-1,7,8,11, 7, 4);

versions = {'parms_v4', 'parms_v2d', 'parms_v2d', 'parms_v2d'};

for p = 1:7
    
    
    for ii = 1:size(mcodes,1)
        
        [output_mainfolder, filename, ~, ~] = get_folder_and_model(mcodes(ii,:));
        
            cd(['C:\Users\u0167448\Documents\GitHub\biophysical-muscle-model\Model output\RMSD\', versions{ii}])
        load([filename, '_RMSD.mat'], 'RMSD', 'AMPs', 'ISIs')
        
        sAMPs = AMPs;
        sISIs = ISIs;
        
        %% remove some trials
        RMSD(1,:,sAMPs < .03,2,:) = nan;
        RMSD(1,sISIs > 1,:,2,:) = nan;
        
        %% average over identified force levels
        sISI = .001;
        sAMP = .0383;
        
            cd(['C:\Users\u0167448\Documents\GitHub\biophysical-muscle-model\Model output\SRS\', versions{ii}])
        load([filename, '_SRS.mat'],'F0', 'AMPs', 'ISIs')
        
%         sACTis,ISIs==sISI, AMPs==sAMP
        
        % average
        for k = iFs
            for i = 1:length(th)-1
                id = F0(:,ISIs == sISI, AMPs == sAMP,k) > th(i) & F0(:,ISIs == sISI, AMPs == sAMP,k) <= th(i+1);
                
                RMSDs(i,:,:,k, p, ii) = mean(RMSD(id,:,:,k,p),1,'omitnan');
                F0s(i,:,:,k, ii) = mean(F0(id,ismember(ISIs, sISIs),ismember(AMPs, sAMPs),k), 1,'omitnan');
            end
        end
    end
    
    %% correct for individual offsets
    % correct for individual offset
    for ii = 1:size(mcodes,1)
        for k = iFs
            RMSDc(:,:,:,k,p,ii) = RMSDs(:,:,:,k,p,ii) - mean(RMSDs(:,:,:,k,p,ii),'all','omitnan') + mean(RMSDs(:,:,:,:,p,ii),'all','omitnan');
        end
    end
    
    %%
end


%% make figure
ps = [7 7 7;
      2 2 2;
      5 5 5];

close all
clc

% yrange = [-.05 .15];

color = get(gca,'colororder');
pcolors = flip(parula(7));
color = [color(2,:); pcolors(4:end-1,:);pcolors(4:end-1,:)];


%%
figure(1)

ms = [7 6 6 5];
sym = {'s','h','v','o'};

for p = 1:3
    % figure;
    for kk = 1:size(mcodes,1)
%         color = get(gca,'colororder');
        
        
        % effect of activation
        AMPid = 7;
        ISIid = 3;
        pCaid = 1:4;
        
        subplot(3,3,1 + (p-1)*3);
         plot(mean(F0s(pCaid,ISIid,AMPid,:, kk), 4, 'omitnan'), mean(RMSDc(pCaid,ISIid,AMPid,:,ps(p,1), kk),4, 'omitnan'), '--', 'color', color(kk,:)); hold on
     
        if kk == 1
            if (p-1)*3 == 0
                yrange = [-.02 .05];
            else
                yrange = [-.02 .2];
            end
            
            plot(ones(1,2) * mean(F0s(pCaid(3),ISIid,AMPid,:, kk), 4, 'omitnan'), yrange, 'k--'); hold on
            
            for ii = pCaid
                plot(ones(1,2) * mean(F0s(ii,ISIid,AMPid,:, kk), 4, 'omitnan'), yrange,'k:')
            end
        end
        
        pCaid = 2:3;
        errorbar(mean(F0s(pCaid,ISIid,AMPid,:, kk), 4, 'omitnan'), mean(RMSDc(pCaid,ISIid,AMPid,:,ps(p,1), kk),4, 'omitnan'), ...
            std(RMSDc(pCaid,ISIid,AMPid,:,ps(p,1),kk),1,4, 'omitnan'), std(RMSDc(pCaid,ISIid,AMPid,:,ps(p,1),kk),1,4, 'omitnan'),...
            std(F0s(pCaid,ISIid,AMPid,:,kk),1,4, 'omitnan'), std(F0s(pCaid,ISIid,AMPid,:,kk),1,4, 'omitnan'),sym{kk}, 'color', color(kk,:),'markerfacecolor', color(kk,:), 'markersize',  ms(kk), 'Capsize',3); hold on
       
        pCaid = [1 4];
        errorbar(mean(F0s(pCaid,ISIid,AMPid,:, kk), 4, 'omitnan'), mean(RMSDc(pCaid,ISIid,AMPid,:,ps(p,1), kk),4, 'omitnan'), ...
            std(RMSDc(pCaid,ISIid,AMPid,:,ps(p,1),kk),1,4, 'omitnan'), std(RMSDc(pCaid,ISIid,AMPid,:,ps(p,1),kk),1,4, 'omitnan'),...
            std(F0s(pCaid,ISIid,AMPid,:,kk),1,4, 'omitnan'), std(F0s(pCaid,ISIid,AMPid,:,kk),1,4, 'omitnan'),sym{kk}, 'color', color(kk,:),'markerfacecolor', [1 1 1], 'markersize', ms(kk), 'Capsize',3); hold on
        
        % effect of amplitude
        AMPid = [2, 3, 4, 7];
        pCaid = 3;
        
        subplot(3,3,2 + (p-1)*3);
           plot(sAMPs(AMPid), squeeze(mean(RMSDc(pCaid,ISIid,AMPid,:,ps(p,2),kk),4, 'omitnan')), '--','color', color(kk,:)); hold on
         
        if kk == 1
            plot(ones(1,2) * 0.0383, yrange, 'k-.'); hold on
        end
        
        AMPid = [2, 3, 4];
        errorbar(sAMPs(AMPid), squeeze(mean(RMSDc(pCaid,ISIid,AMPid,:,ps(p,2),kk),4, 'omitnan')), ...
            squeeze(std(RMSDc(pCaid,ISIid,AMPid,:,ps(p,2),kk),1,4, 'omitnan')),sym{kk}, 'color', color(kk,:),'markerfacecolor', [1 1 1], 'markersize', ms(kk), 'Capsize',3); hold on
        
        AMPid = 7;
        errorbar(sAMPs(AMPid), squeeze(mean(RMSDc(pCaid,ISIid,AMPid,:,ps(p,2),kk),4, 'omitnan')), ...
            squeeze(std(RMSDc(pCaid,ISIid,AMPid,:,ps(p,2),kk),1,4, 'omitnan')),sym{kk}, 'color', color(kk,:),'markerfacecolor', color(kk,:), 'markersize', ms(kk), 'Capsize',3); hold on
        
        xlim([0 .0383])
        
        % effect of ISI
        ISIid = [1, 3:5, 7];
        AMPid = 7;
        
        subplot(3,3,3 + (p-1)*3);
              plot(sISIs(ISIid), squeeze(mean(RMSDc(pCaid,ISIid,AMPid,:,ps(p,3),kk),4, 'omitnan')), '--','color', color(kk,:)); hold on
          
        if kk == 1
            plot(ones(1,2) * .1, yrange, 'k-.'); hold on
        end
        
        ISIid = [5,7];
        errorbar(sISIs(ISIid), squeeze(mean(RMSDc(pCaid,ISIid,AMPid,:,ps(p,3),kk),4, 'omitnan')), ...
            squeeze(std(RMSDc(pCaid,ISIid,AMPid,:,ps(p,3),kk),1,4, 'omitnan')),sym{kk}, 'color', color(kk,:),'markerfacecolor', [1 1 1], 'markersize',ms(kk), 'Capsize',3); hold on
        
        ISIid = [1, 3, 4];
        errorbar(sISIs(ISIid), squeeze(mean(RMSDc(pCaid,ISIid,AMPid,:,ps(p,3),kk),4, 'omitnan')), ...
            squeeze(std(RMSDc(pCaid,ISIid,AMPid,:,ps(p,3),kk),1,4, 'omitnan')),sym{kk}, 'color', color(kk,:),'markerfacecolor', color(kk,:), 'markersize', ms(kk), 'Capsize',3); hold on
        
        set(gca, 'XScale', 'log', 'Xlim', [1e-3 1e1])
        
        
    end
end


% make nice
figure(1)

for j = 1:9
    subplot(3,3,j)
      set(gca, 'Fontsize', 6)
    
    box off
    %     title(titles{j})
    %     xlabel(xlabels{j})
    %     ylabel(ylabels)
    ylim(yrange)

    if j < 7
        set(gca,'Xticklabel', [])
    end
    
    if j ~= 1 && j ~= 4 && j ~= 7
        set(gca,'Yticklabel', [])
        set(gca, 'Ycolor', 'none')
    else
        ylabel('RMSD (F_{0})',  'Fontsize', 8)
        xlim([0 1.02])
    end
    
    if j == 7
        %         title('Overall')
        xlabel('Activation (F_{0})',  'Fontsize', 8)
    elseif j == 8
        xlabel('Amplitude (L_0)',  'Fontsize', 8)
    elseif j == 9
        xlabel('Recovery time (s)',  'Fontsize', 8)
    end
    
%     if 
    
end

%% titles
subtitles = {'Overall', 'Pre-stretch', 'Test stretch'};

figure(1)
subplot(331)
title('Effect of activation', 'fontsize', 8)
subtitle(subtitles{1}, 'fontsize', 8)

subplot(332)
subtitle(subtitles{1}, 'fontsize', 8)
title('Effect of amplitude', 'fontsize', 8)

subplot(333)
subtitle(subtitles{1}, 'fontsize', 8)
title('Effect of recovery', 'fontsize', 8)

subplot(334)
% subtitle('Effect of activation', 'fontsize', 8)
subtitle(subtitles{2}, 'fontsize', 8)

subplot(335)
% subtitle('Effect of amplitude', 'fontsize', 8)
subtitle(subtitles{2}, 'fontsize', 8)

subplot(336)
% subtitle('Effect of recovery', 'fontsize', 8)
subtitle(subtitles{2}, 'fontsize', 8)

subplot(337)
% subtitle('Effect of activation', 'fontsize', 8)
subtitle(subtitles{3}, 'fontsize', 8)

subplot(338)
% subtitle('Effect of amplitude', 'fontsize', 8)
subtitle(subtitles{3}, 'fontsize', 8)

subplot(339)
% subtitle('Effect of recovery', 'fontsize', 8)
subtitle(subtitles{3}, 'fontsize', 8)

%% make manual legend
% if ishandle(2), close(2); end; figure(2)
modelnames = {'Hill model','XB model','XB coop','XB coop + FD'};
figure(1)
ys = .19;

subplot(331)
for i = 1:4
    fake_data = -i * .02  + ys * [.95 1 1.05];
    errorbar(.1, mean(fake_data, 2), std(fake_data,1,2),sym{i}, 'color', color(i,:),'markerfacecolor', color(i,:),'markersize', ms(i)*.5, 'Capsize',3); hold on
    errorbar(.8, mean(fake_data, 2), std(fake_data,1,2),sym{i}, 'color', color(i,:), 'markerfacecolor', [1 1 1],'markersize', ms(i)*.5, 'Capsize',3); hold on

    text(.45, ys - i * .02, modelnames{i}, 'fontsize', 7, 'horizontalalignment', 'center')

end

%%
text(.1, ys, 'HD', 'fontsize', 7, 'horizontalalignment', 'center')
text(.8, ys, 'HD', 'fontsize', 7, 'horizontalalignment', 'center')

plot(.8, ys, 'rx', 'linewidth', 2)

%%
figure(1)
set(gcf,'units','centimeters','position',[10 10 19 11])


%% A, B labels
subplot(331)
text(-.1, yrange(2)*1.15, 'A', 'fontsize', 12,'fontweight', 'bold')
subplot(334)
text(-.1, yrange(2)*1.15, 'B', 'fontsize', 12,'fontweight', 'bold')
subplot(337)
text(-.1, yrange(2)*1.15, 'C', 'fontsize', 12,'fontweight', 'bold')

%% optionally export to PNG

if savefig
cd(['C:\Users\',username,'\OneDrive\9. Short-range stiffness\figures\MAT'])
       
figure(1)
exportgraphics(gcf,['Fig7.png'])
end
