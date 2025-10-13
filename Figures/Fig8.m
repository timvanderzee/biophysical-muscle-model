clear all; close all; clc
savefig = 1;

[username, githubfolder] = get_paths();


%% Model
tid = [3 1 7];


showbar = 0;
showline = 1;
% figure(1)
filenames = {'Hill_regular_SRS', 'biophysical_no_regular_SRS', 'biophysical_full_regular_SRS', 'biophysical_full_alternative_SRS'};
% visualize_model_SRS_simple(filenames, th, id)

th = [0 .07 .15 .3 .5 .7 1.5];
tid(1) = 4;

iFs = [2,3,5,6,7,8,11];
yrange = [.3 1.5];

figure(1)
color = get(gca,'colororder');
pcolors = flip(parula(7));
color = [color(2,:); pcolors(4:end-1,:);pcolors(4:end-1,:)];

for kk = 1:length(filenames)
    
    load(filenames{kk},'Stest', 'Scond', 'AMPs', 'pCas', 'ISIs', 'F0')

    % average
    SRSrel_m = nan(length(th)-1,length(ISIs),length(AMPs),iFs(end));
    F0s_m = nan(length(th)-1, length(ISIs),length(AMPs),iFs(end));
    
    for k = iFs
        for i = 1:length(th)-1
            id = F0(:,1,7,k) > th(i) & F0(:,1,7,k) <= th(i+1);
            
            SRSrel_m(i,:,:,k) = mean(Stest(id,:,:,k),1,'omitnan') ./ mean(Scond(id,1,7,k),'all', 'omitnan');
            F0s_m(i,:,:,k) = mean(F0(id,:,:,k), 1,'omitnan');
        end
    end
    
    figure(1)
    if kk == 1
        subplot(131)
        plot(ones(1,2) * mean(F0s_m(tid(1),tid(2),tid(3),:), 4, 'omitnan'), yrange, 'k--'); hold on
        
        subplot(132)
         plot(ones(1,2) * AMPs(tid(3)), yrange, 'k--'); hold on
        
         subplot(133)
         plot(ones(1,2) * ISIs(tid(2)), yrange, 'k--'); hold on
         
    end
    
    subplot(131);
    plot(mean(F0s_m(:,tid(2),tid(3),:), 4, 'omitnan'), mean(SRSrel_m(1:end,tid(2),tid(3),:),4, 'omitnan'), 'color', color(kk,:), 'linewidth',2); hold on
    
    subplot(132);   
    plot(AMPs, squeeze(mean(SRSrel_m(tid(1),tid(2),1:end,:),4, 'omitnan')),'color',color(kk,:), 'linewidth', 2); hold on
    
    subplot(133);   
    plot(ISIs, squeeze(mean(SRSrel_m(tid(1),1:end,tid(3),:),4, 'omitnan')),'color',color(kk,:),'linewidth',2); hold on
   
    set(gca, 'XScale', 'log', 'Xlim', [5e-4 2e1])

   
end

figure(1)

% make nice
titles = {'Effect of activation', 'Effect of amplitude','Effect of recovery'};
xlabels = {'Activation (F_0)', 'Amplitude (L_0)', 'Recovery time (s)'};
ylabels = 'Relative short-range stiffness';

for j = 1:3
    subplot(1,3,j)
    set(gca,'Fontsize', 6)
    box off
    title(titles{j},'Fontsize', 8)
    xlabel(xlabels{j},'Fontsize', 8)
    
    if j == 1
        ylabel(ylabels,'Fontsize', 8)
    end
    
    ylim(yrange)
        yline(1, 'k-')
    
    %     for ii = 1:length(AMPs)
    %         yline(mean(SRSrel(pCaid,ISIid,ii,:),4, 'omitnan'), ':','color', color(ii,:))
    %     end
end


subplot(131)
% xline(mean(F0s(pCaid,jj,ii,:),4, 'omitnan'), 'k:')

subplot(132)
xlim([-.005 .06])
% xline(AMPs(AMPid), 'k:')

subplot(133)
% xline(ISIs(ISIid), 'k:')

xlim([5e-4 2e1])


%% Data
eth = [0 .07 .25 .7 1.5];
tid = [3 1 7];

Cid = tid(1);
ISIid = tid(2);
AMPid = tid(3);

% [username, githubfolder] = get_paths();
% cd([githubfolder, '\biophysical-muscle-model\Data'])
load('SRS_data_v3.mat', 'SRS_pre', 'F0', 'SRS_post')

% average some pCas
iFs = 1:11;
%  iFs = [1 2 3, 5, 6, 7, 8, 10, 11];
SRSrel = nan(length(eth)-1,7,8,length(iFs));
F0s = nan(length(eth)-1, 7,8,length(iFs));

for k = iFs
    for i = 1:length(eth)-1
        id = F0(:,1,7,k) > eth(i) & F0(:,1,7,k) <= eth(i+1);
        
        SRSrel(i,:,:,k) = mean(SRS_post(id,:,:,k),1,'omitnan') ./ mean(SRS_pre(id,:,7,k),'all', 'omitnan');
%         SRSrel(i,:,:,k) = mean(SRS_post(id,:,:,k),1,'omitnan') ./ mean(SRS_post(id,1,1,k),'all', 'omitnan');
        F0s(i,:,:,k) = mean(F0(id,:,:,k), 1,'omitnan');
    end
end


% summary plot
eAMPs = [0 12 38 121 216 288 383 682]/10000;
eISIs = [1 10 100 316 1000 3160 10000]/1000;

iid = [1 3 4 5 7];
aid = [1, 3, 4, 7];
pCaid = 2:3;

subplot(131)
SST(1) = sum((squeeze(mean(SRSrel(:,ISIid, AMPid, :),4,'omitnan')) - mean(mean(SRSrel(:,ISIid, AMPid, :),4,'omitnan'),'omitnan')).^2,'omitnan');

plot(squeeze(F0s(:,1,7,:)), squeeze(SRSrel(:,ISIid,AMPid,:)), '.', 'color', [.5 .5 .5]); hold on

errorbar(mean(F0s(pCaid,1,7,:),4,'omitnan'), mean(SRSrel(pCaid,ISIid,AMPid,:),4,'omitnan'),...
    std(SRSrel(pCaid,ISIid,AMPid,:),1,4,'omitnan'),std(SRSrel(pCaid,ISIid,AMPid,:),1,4,'omitnan'),...
    std(F0s(pCaid,1,7,:),1,4,'omitnan'),std(F0s(pCaid,1,7,:),1,4,'omitnan'),'o', 'color', [.5 .5 .5], 'markerfacecolor', [.5 .5 .5]); hold on

pCaid = [1 4];
errorbar(mean(F0s(pCaid,1,7,:),4,'omitnan'), mean(SRSrel(pCaid,ISIid,AMPid,:),4,'omitnan'),...
    std(SRSrel(pCaid,ISIid,AMPid,:),1,4,'omitnan'),std(SRSrel(pCaid,ISIid,AMPid,:),1,4,'omitnan'),...
    std(F0s(pCaid,1,7,:),1,4,'omitnan'),std(F0s(pCaid,1,7,:),1,4,'omitnan'),'o', 'color', [.5 .5 .5], 'markerfacecolor', [1 1 1]); hold on

subplot(132)
SST(2) = sum((squeeze(mean(SRSrel(Cid,ISIid, aid, :),4,'omitnan')) - mean(mean(SRSrel(Cid,ISIid, aid, :),4,'omitnan'),'omitnan')).^2,'omitnan');

plot(eAMPs(aid), squeeze(SRSrel(Cid,ISIid,aid,:)), '.', 'color', [.5 .5 .5]); hold on
errorbar(eAMPs(aid(1:3)), squeeze(mean(SRSrel(Cid,ISIid,aid(1:3),:), 4, 'omitnan')), squeeze(std(SRSrel(Cid,ISIid,aid(1:3),:), 1, 4, 'omitnan')), 'o', 'color', [.5 .5 .5], 'markerfacecolor', [1 1 1])
errorbar(eAMPs(aid(4)), squeeze(mean(SRSrel(Cid,ISIid,aid(4),:), 4, 'omitnan')), squeeze(std(SRSrel(Cid,ISIid,aid(4),:), 1, 4, 'omitnan')), 'o', 'color', [.5 .5 .5], 'markerfacecolor', [.5 .5 .5])

subplot(133)
SST(3) = sum((squeeze(mean(SRSrel(Cid,iid, AMPid, :),4,'omitnan')) - mean(mean(SRSrel(Cid,iid, AMPid, :),4,'omitnan'),'omitnan')).^2,'omitnan');

semilogx(eISIs(iid), squeeze(SRSrel(Cid,iid,AMPid,:)), '.', 'color', [.5 .5 .5]); hold on
errorbar(eISIs(iid(4:5)), squeeze(mean(SRSrel(Cid,iid(4:5),AMPid,:), 4, 'omitnan')), squeeze(std(SRSrel(Cid,iid(4:5),AMPid,:), 1, 4, 'omitnan')), 'o', 'color', [.5 .5 .5], 'markerfacecolor', [1 1 1])
errorbar(eISIs(iid(1:3)), squeeze(mean(SRSrel(Cid,iid(1:3),AMPid,:), 4, 'omitnan')), squeeze(std(SRSrel(Cid,iid(1:3),AMPid,:), 1, 4, 'omitnan')), 'o', 'color', [.5 .5 .5], 'markerfacecolor', [.5 .5 .5])

%%
figure(1)
set(gcf,'units','normalized','position',[.1 .1 .45 .25])


%% A, B labels
figure(1)
subplot(131)
text(-.15, 1.55, 'A', 'fontsize', 12,'fontweight','bold')

subplot(132)
text(-.008, 1.55, 'B', 'fontsize', 12,'fontweight','bold')

subplot(133)
text(1e-4, 1.55, 'C', 'fontsize', 12,'fontweight','bold')
    

%% calc and diplay R2
% selected conditions
sAMPs = [0, 0.0012, 0.0038, 0.0121, 0.0383];
sISIs = [0.001, 0.1, 0.3160, 1, 10];
sACTis = [1 2 4 6];

sACTi = 3;
sISI = .001;
sAMP = .0383;

for kk = 1:length(filenames)
    
    load(filenames{kk},'Stest', 'Scond', 'AMPs', 'pCas', 'ISIs', 'F0')

    % average
    SRSrel_m = nan(length(th)-1,length(ISIs),length(AMPs),iFs(end));
    F0s_m = nan(length(th)-1, length(ISIs),length(AMPs),iFs(end));
    
    for k = iFs
        for i = 1:length(th)-1
            id = F0(:,1,7,k) > th(i) & F0(:,1,7,k) <= th(i+1);
            
            SRSrel_m(i,:,:,k) = mean(Stest(id,:,:,k),1,'omitnan') ./ mean(Scond(id,1,7,k),'all', 'omitnan');
            F0s_m(i,:,:,k) = mean(F0(id,:,:,k), 1,'omitnan');
        end
    end
    
    % activation
    SSE(kk,1) = sum((squeeze(mean(SRSrel(:,eISIs==sISI, eAMPs==sAMP, :),4,'omitnan')) - squeeze(mean(SRSrel_m(sACTis,ISIs==sISI, AMPs==sAMP, :),4,'omitnan'))).^2,'omitnan');
    
    % amplitude
    SSE(kk,2) = sum((squeeze(mean(SRSrel(sACTi,eISIs==sISI, ismember(eAMPs, sAMPs), :),4,'omitnan')) - squeeze(mean(SRSrel_m(sACTis(sACTi),ISIs==sISI, ismember(AMPs, sAMPs), :),4,'omitnan'))).^2,'omitnan');
    
    % recovery time
    SSE(kk,3) = sum((squeeze(mean(SRSrel(sACTi,ismember(eISIs, sISIs), eAMPs==sAMP, :),4,'omitnan')) - squeeze(mean(SRSrel_m(sACTis(sACTi),ismember(ISIs, sISIs), AMPs==sAMP, :),4,'omitnan'))).^2,'omitnan');
end

R2 = 1 - SSE ./ SST;

%% legend
figure(1)
subplot(131)

xl = [.42 .5; .001 .005; 1.5e-3 4e-3];
yl = [1.35; 0.6; 1.35];

r = [5 5 100];
dy = .07;

for i = 1:3
    subplot(1,3,i)
    
    plot(xl(i,:), [yl(i) yl(i)], 'color', color(1,:),'linewidth',2)
    plot(xl(i,:), [yl(i) yl(i)]-dy, 'color', color(2,:),'linewidth',2)
    plot(xl(i,:), [yl(i) yl(i)]-2*dy, 'color', color(3,:),'linewidth',2)
    plot(xl(i,:), [yl(i) yl(i)]-3*dy, 'color', color(4,:),'linewidth',2)
    

    text(xl(i,1), yl(i)+dy, 'Model', 'fontsize',8)
    text(xl(i,2) + diff(xl(i,:)) * r(i), yl(i)+dy,'R^2', 'fontsize',8,'horizontalalignment','center')
    
    text(xl(i,2) + diff(xl(i,:)) * .1, yl(i), 'Hill', 'fontsize',6)
    text(xl(i,2) + diff(xl(i,:)) * r(i), yl(i), num2str(round(R2(1,i),2),3), 'fontsize',6,'horizontalalignment','center')

    text(xl(i,2) + diff(xl(i,:)) * .1, yl(i)-dy, 'XB', 'fontsize',6)
    text(xl(i,2) + diff(xl(i,:)) * r(i), yl(i)-dy, num2str(round(R2(2,i),2),3), 'fontsize',6,'horizontalalignment','center')

    text(xl(i,2) + diff(xl(i,:)) * .1, yl(i)-2*dy, 'XB coop', 'fontsize',6)
    text(xl(i,2) + diff(xl(i,:)) * r(i), yl(i)-2*dy, num2str(round(R2(3,i),2),3), 'fontsize',6,'horizontalalignment','center')

    text(xl(i,2) + diff(xl(i,:)) * .1, yl(i)-3*dy, 'XB coop + FD', 'fontsize',6)
    text(xl(i,2) + diff(xl(i,:)) * r(i), yl(i)-3*dy, num2str(round(R2(4,i),2),3), 'fontsize',6,'horizontalalignment','center')
end

subplot(131)
text(.9, 0.5, 'HD', 'horizontalalignment','center', 'fontsize',8)
text(1, 0.5, 'HD', 'horizontalalignment','center', 'fontsize',8)
plot(1, 0.5, 'rx','linewidth',2)

errorbar(.9, .4, .05, 'o', 'color', [.5 .5 .5], 'markersize', 5, 'markerfacecolor', [.5 .5 .5])
errorbar(1, .4, .05, 'o', 'color', [.5 .5 .5], 'markersize', 5, 'markerfacecolor', [1 1 1])

xlim([0 1.05])

subplot(132)
xlim([-.001 .05])

%% save and export

if savefig
cd(['C:\Users\',username,'\OneDrive\9. Short-range stiffness\figures\MAT'])
       
figure(1)
exportgraphics(gcf,'Fig8.png')
end

return

%% calc AIC

    % include all trials
%     for iii = 1:3
%         
%         if iii == 3 % all trials
%             iid = [1 3 4 5 7];
%             aid = [1 2 3 4 7];
%             y = SRSrel;
%             
%         elseif iii == 2 % no HD
%             iid = [1 3 4 5 7];
%             aid = [1 2 3 4 7];
%             
%             % exclude high history dependent trials
%             y = SRSrel;
%             y(:,[1 3 4], 7,:) = nan;
%             
%         elseif iii == 1 % HD
%             y = SRSrel;
% 
%             iid = [1 3 4];
%             aid = 7; 
%         end
%         
%         % over all conditions in the dataset
%         mSRSrel = mean(y, 4, 'omitnan');
%         N(kk,iii) = sum(isfinite(mSRSrel(:,iid,aid)),'all','omitnan');
%         mSRSrel_m = mean(SRSrel_m, 4, 'omitnan');
%         mmSRSrel = mean(mSRSrel(:,iid,aid), 'all', 'omitnan');
% 
% %         oSST(1,iii) = sum((mSRSrel(:,iid,aid) - mmSRSrel).^2,'all','omitnan');
% %         oSSE(kk,iii) = sum((mSRSrel(:,iid,aid) - mSRSrel_m(:,iid,aid)).^2,'all','omitnan');
% 
%     end
    
    % for each trial individually
%     oSSEs(:,:,:,kk) = (mSRSrel - mSRSrel_m).^2;
% 
%     
%     atitles = {'Passive', 'Low activation', 'Medium activation', 'Max. activation'};
%     for iii = 1:4
%         figure(2)
%         subplot(4,4,iii + (kk-1)*4)
%         surf(aAMPs(iid,aid)*100, aISIs(iid,aid), squeeze((mSRSrel(iii,iid,aid) - mSRSrel_m(iii,iid,aid)).^2))
%         set(gca,'YScale','log', 'Clim', [0 .1])
%         zlim([0 .2])
%         
%         yt = get(gca,'Ytick');
%         set(gca,'Yticklabel',yt);
%         
%         if kk == 1
%             title(atitles{iii})
%         end
%         
%         if iii > 1
%             set(gca,'Zticklabel', [])
%         end
%         
%         if kk < 4
%             set(gca,'Xticklabel', [], 'Yticklabel', [])
%         end
%     end

%% R2
clc

R2 = 1 - SSE ./ SST;

% n = repmat([4 4 5], 4, 1); % number of data points
% k = repmat([5; 8; 10; 12], 1, 3); % number of parameters

% AIC = 2*k + n.*log(SSE./n);

% AICc = AIC + (2*k.^2 + 2*k) ./ (n-k-1)

%% first sum, then compute AIC
tAIC = 2*k(:,1) + sum(n,2) .* log(sum(SSE,2)./sum(n,2));

%% overal R2 and AIC
id = 3;

R2 = 1 - oSSE(:,id) ./ oSST(id);

% n = numel(mSRSrel_m);
n = N(1,id);

k = [5 8 10 12]';

AIC = 2*k + n.*log(oSSE(:,id)./n);

%% visualize
close all
figure(10)

% set(gca,'Colororder',color)
b = bar(oSSE'./N');

for i = 1:4
    set(b(i), 'FaceColor', color(i,:))
end

set(gca, 'Xticklabel', {'HD trials', 'No HD trials', 'All trials'})
ylabel('Mean SRS error')
box off
legend('Hill','XB no coop', 'XB coop', 'XB coop + FD', 'location', 'best')
legend boxoff




