clear all; close all; clc
[username, githubfolder] = get_paths();
th = [0 .07 .25 .7 1.5];
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
        
        SRSrel(i,:,:,k) = mean(SRS_post(id,:,:,k),1,'omitnan') ./ mean(SRS_pre(id,:,7,k),'all', 'omitnan');
%         SRSrel(i,:,:,k) = mean(SRS_post(id,:,:,k),1,'omitnan') ./ mean(SRS_post(id,1,1,k),'all', 'omitnan');
        F0s(i,:,:,k) = mean(F0(id,:,:,k), 1,'omitnan');
    end
end


% summary plot
AMPs = [0 12 38 121 216 288 383 682]/10000;
ISIs = [1 10 100 316 1000 3160 10000]/1000;

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

plot(AMPs(aid), squeeze(SRSrel(Cid,ISIid,aid,:)), '.', 'color', [.5 .5 .5]); hold on
errorbar(AMPs(aid(1:3)), squeeze(mean(SRSrel(Cid,ISIid,aid(1:3),:), 4, 'omitnan')), squeeze(std(SRSrel(Cid,ISIid,aid(1:3),:), 1, 4, 'omitnan')), 'o', 'color', [.5 .5 .5], 'markerfacecolor', [1 1 1])
errorbar(AMPs(aid(4)), squeeze(mean(SRSrel(Cid,ISIid,aid(4),:), 4, 'omitnan')), squeeze(std(SRSrel(Cid,ISIid,aid(4),:), 1, 4, 'omitnan')), 'o', 'color', [.5 .5 .5], 'markerfacecolor', [.5 .5 .5])


subplot(133)
SST(3) = sum((squeeze(mean(SRSrel(Cid,iid, AMPid, :),4,'omitnan')) - mean(mean(SRSrel(Cid,iid, AMPid, :),4,'omitnan'),'omitnan')).^2,'omitnan');

semilogx(ISIs(iid), squeeze(SRSrel(Cid,iid,AMPid,:)), '.', 'color', [.5 .5 .5]); hold on
errorbar(ISIs(iid(4:5)), squeeze(mean(SRSrel(Cid,iid(4:5),AMPid,:), 4, 'omitnan')), squeeze(std(SRSrel(Cid,iid(4:5),AMPid,:), 1, 4, 'omitnan')), 'o', 'color', [.5 .5 .5], 'markerfacecolor', [1 1 1])
errorbar(ISIs(iid(1:3)), squeeze(mean(SRSrel(Cid,iid(1:3),AMPid,:), 4, 'omitnan')), squeeze(std(SRSrel(Cid,iid(1:3),AMPid,:), 1, 4, 'omitnan')), 'o', 'color', [.5 .5 .5], 'markerfacecolor', [.5 .5 .5])

%% Model
aAMPs = repmat(AMPs, 7, 1);
aISIs = repmat(ISIs(:), 1, 8);
showbar = 0;
showline = 1;
% figure(1)
filenames = {'Hill_regular_SRS', 'biophysical_no_regular_SRS', 'biophysical_full_regular_SRS', 'biophysical_full_alternative_SRS'};
% visualize_model_SRS_simple(filenames, th, id)
for kk = 1:length(filenames)
    
    load(filenames{kk},'Stest', 'Scond', 'AMPs', 'iFs', 'pCas', 'ISIs', 'F0')
    % SRSrel = Stest./Scond(:,:,AMPs == .0383,:);
    % F0s = F0;
    
    iid = [1 3 4 5 7];
    aid = [1, 3, 4, 7];
    
    % average
%     iFs = [1 2 3, 5, 6, 7, 8, 10, 11];
    iFs = [2,3,5,6,7,8,11];
    SRSrel_m = nan(length(th)-1,7,8,iFs(end));
    F0s_m = nan(length(th)-1, 7,8,iFs(end));
    
    for k = iFs
        for i = 1:length(th)-1
            id = F0(:,1,7,k) > th(i) & F0(:,1,7,k) <= th(i+1);
            
            SRSrel_m(i,:,:,k) = mean(Stest(id,:,:,k),1,'omitnan') ./ mean(Scond(id,1,7,k),'all', 'omitnan');
            F0s_m(i,:,:,k) = mean(F0(id,:,:,k), 1,'omitnan');
        end
    end
    
    % plot for all fibers
    yrange = [.3 1.5];
    
    color = get(gca,'colororder');
    pcolors = flip(parula(7));
    color = [color(2,:); pcolors(4:end-1,:);pcolors(4:end-1,:)];

    ms = [7 6 6 5];
    sym = 'shvo';
    
    % effect of activation
    AMPid = tid(3);
    ISIid = tid(2);
    pCaid = 1:4;

    SSE(kk,1) = sum((squeeze(mean(SRSrel(pCaid,ISIid, AMPid, :),4,'omitnan')) - squeeze(mean(SRSrel_m(pCaid,ISIid, AMPid, :),4,'omitnan'))).^2,'omitnan');
    
    figure(1)
    subplot(131);
    
    if kk == 1
        plot(ones(1,2) * mean(F0s_m(pCaid(tid(1)),ISIid,AMPid,:), 4, 'omitnan'), yrange, 'k--'); hold on
        
%         for ii = pCaid
%             plot(ones(1,2) * mean(F0s_m(ii,ISIid,AMPid,:), 4, 'omitnan'), yrange,'k:')
%         end
    end
     
    
    if showbar
    pCaid = 2:3;
    errorbar(mean(F0s_m(pCaid,ISIid,AMPid,:), 4, 'omitnan'), mean(SRSrel_m(pCaid,ISIid,AMPid,:),4, 'omitnan'), ...
        std(SRSrel_m(pCaid,ISIid,AMPid,:),1,4, 'omitnan'), std(SRSrel_m(pCaid,ISIid,AMPid,:),1,4, 'omitnan'),...
        std(F0s_m(pCaid,ISIid,AMPid,:),1,4, 'omitnan'), std(F0s_m(pCaid,ISIid,AMPid,:),1,4, 'omitnan'),sym(kk), 'color', color(kk,:),'markerfacecolor', color(kk,:), 'markersize', ms(kk)); hold on
    
    pCaid = [1 4];
    errorbar(mean(F0s_m(pCaid,ISIid,AMPid,:), 4, 'omitnan'), mean(SRSrel_m(pCaid,ISIid,AMPid,:),4, 'omitnan'), ...
        std(SRSrel_m(pCaid,ISIid,AMPid,:),1,4, 'omitnan'), std(SRSrel_m(pCaid,ISIid,AMPid,:),1,4, 'omitnan'),...
        std(F0s_m(pCaid,ISIid,AMPid,:),1,4, 'omitnan'), std(F0s_m(pCaid,ISIid,AMPid,:),1,4, 'omitnan'),sym(kk), 'color', color(kk,:),'markerfacecolor', [1 1 1], 'markersize', ms(kk)); hold on
    end
    
    if showline
        pCaid = 1:4;
        plot(mean(F0s_m(pCaid,ISIid,AMPid,:), 4, 'omitnan'), mean(SRSrel_m(pCaid,ISIid,AMPid,:),4, 'omitnan'), 'color', color(kk,:), 'linewidth',2); hold on
    end
    
    % effect of amplitude
    AMPid = aid(1:3);
    pCaid = tid(1);
    
    subplot(132);
    SSE(kk,2) = sum((squeeze(mean(SRSrel(pCaid,ISIid, aid, :),4,'omitnan')) - squeeze(mean(SRSrel_m(pCaid,ISIid, aid, :),4,'omitnan'))).^2,'omitnan');
    
    if kk == 1
        plot(ones(1,2) * AMPs(tid(3)), yrange, 'k--'); hold on
    end
    
    if showbar
    errorbar(AMPs(AMPid), squeeze(mean(SRSrel_m(pCaid,ISIid,AMPid,:),4, 'omitnan')), ...
        squeeze(std(SRSrel_m(pCaid,ISIid,AMPid,:),1,4, 'omitnan')),sym(kk), 'color', color(kk,:),'markerfacecolor', [1 1 1], 'markersize', ms(kk)); hold on
    
    AMPid = [aid(4) 8]; % add the largest
    errorbar(AMPs(AMPid), squeeze(mean(SRSrel_m(pCaid,ISIid,AMPid,:),4, 'omitnan')), ...
        squeeze(std(SRSrel_m(pCaid,ISIid,AMPid,:),1,4, 'omitnan')),sym(kk), 'color', color(kk,:),'markerfacecolor', color(kk,:), 'markersize', ms(kk)); hold on
    end
    if showline
        AMPid = 1:length(AMPs);
        plot(AMPs(AMPid), squeeze(mean(SRSrel_m(pCaid,ISIid,AMPid,:),4, 'omitnan')),'color',color(kk,:), 'linewidth', 2); hold on
    end
    
    % effect of ISI
    ISIid = iid(4:5);
    AMPid = tid(3);
    
    subplot(133);
    SSE(kk,3) = sum((squeeze(mean(SRSrel(pCaid,iid, AMPid, :),4,'omitnan')) - squeeze(mean(SRSrel_m(pCaid,iid, AMPid, :),4,'omitnan'))).^2,'omitnan');
        
    if kk == 1
        plot(ones(1,2) * ISIs(tid(2)), yrange, 'k--'); hold on
    end
    
    if showbar
    errorbar(ISIs(ISIid), squeeze(mean(SRSrel_m(pCaid,ISIid,AMPid,:),4, 'omitnan')), ...
        squeeze(std(SRSrel_m(pCaid,ISIid,AMPid,:),1,4, 'omitnan')),sym(kk), 'color', color(kk,:),'markerfacecolor', [1 1 1], 'markersize',ms(kk)); hold on
    
    ISIid = iid(1:3);
    errorbar(ISIs(ISIid), squeeze(mean(SRSrel_m(pCaid,ISIid,AMPid,:),4, 'omitnan')), ...
        squeeze(std(SRSrel_m(pCaid,ISIid,AMPid,:),1,4, 'omitnan')),sym(kk), 'color', color(kk,:),'markerfacecolor', color(kk,:), 'markersize', ms(kk)); hold on
    
    end
    if showline
        ISIid = 1:length(ISIs);
        plot(ISIs(ISIid), squeeze(mean(SRSrel_m(pCaid,ISIid,AMPid,:),4, 'omitnan')),'color',color(kk,:),'linewidth',2); hold on
    end
    
    set(gca, 'XScale', 'log', 'Xlim', [5e-4 2e1])

    % include all trials
    for iii = 1:3
        
        if iii == 3 % all trials
            iid = [1 3 4 5 7];
            aid = [1 2 3 4 7];
            y = SRSrel;
            
        elseif iii == 2 % no HD
            iid = [1 3 4 5 7];
            aid = [1 2 3 4 7];
            
            % exclude high history dependent trials
            y = SRSrel;
            y(:,[1 3 4], 7,:) = nan;
            
        elseif iii == 1 % HD
            y = SRSrel;

            iid = [1 3 4];
            aid = 7; 
        end
        
        % over all conditions in the dataset
        mSRSrel = mean(y, 4, 'omitnan');
        N(kk,iii) = sum(isfinite(mSRSrel(:,iid,aid)),'all','omitnan');
        mSRSrel_m = mean(SRSrel_m, 4, 'omitnan');
        mmSRSrel = mean(mSRSrel(:,iid,aid), 'all', 'omitnan');

        oSST(1,iii) = sum((mSRSrel(:,iid,aid) - mmSRSrel).^2,'all','omitnan');
        oSSE(kk,iii) = sum((mSRSrel(:,iid,aid) - mSRSrel_m(:,iid,aid)).^2,'all','omitnan');

    end
    
    % for each trial individually
    oSSEs(:,:,:,kk) = (mSRSrel - mSRSrel_m).^2;

    
    atitles = {'Passive', 'Low activation', 'Medium activation', 'Max. activation'};
    for iii = 1:4
        figure(2)
        subplot(4,4,iii + (kk-1)*4)
        surf(aAMPs(iid,aid)*100, aISIs(iid,aid), squeeze((mSRSrel(iii,iid,aid) - mSRSrel_m(iii,iid,aid)).^2))
        set(gca,'YScale','log', 'Clim', [0 .1])
        zlim([0 .2])
        
        yt = get(gca,'Ytick');
        set(gca,'Yticklabel',yt);
        
        if kk == 1
            title(atitles{iii})
        end
        
        if iii > 1
            set(gca,'Zticklabel', [])
        end
        
        if kk < 4
            set(gca,'Xticklabel', [], 'Yticklabel', [])
        end
    end
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

%%
figure(1)
set(gcf,'units','normalized','position',[.1 .1 .45 .3])

%%
cd(['C:\Users\',username,'\OneDrive\9. Short-range stiffness\figures\MAT'])
       
figure(1)
exportgraphics(gcf,['Fig8.png'])

return
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




