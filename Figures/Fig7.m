clear all; close all; clc
[username, githubfolder] = get_paths();
mcodes = [2 1 1; 1 1 3; 1 1 1; 1 2 1];

iFs = [1 2 3, 5, 6, 7, 8, 10, 11];

% iFs = [2,3,5,6,7,8,11];

th = [0 .07 .25 .7 1.5];

RMSDs = nan(length(th)-1,7,8,length(iFs), 7, 2);
F0s = nan(length(th)-1, 7,8,length(iFs), 2);
RMSDc = nan(length(th)-1,7,8,length(iFs), 7, 2);
% ps = [1, 2, 7;
%       4, 5, ];

% ps = [4, 5, 7];

for p = 1:7
    
    
    for ii = 1:size(mcodes,1)
        [output_mainfolder, filename, ~, ~] = get_folder_and_model(mcodes(ii,:));
        
        load([filename, '_RMSD.mat'], 'RMSD', 'AMPs', 'ISIs')
        
        %% remove some trials
        RMSD(1,:,AMPs < .03,2,:) = nan;
        RMSD(1,ISIs > 1,:,2,:) = nan;
        
        % %% one fiber, all conditions
        % close all
        % aAMPs = repmat(AMPs, 7, 1);
        % aISIs = repmat(ISIs(:), 1, 8);
        %
        % iF = 2;
        % for i = 1:size(RMSD,1)
        %     eps = squeeze(RMSD(i,:,:,iF));
        %
        %     figure(1)
        %     nexttile
        %     surf(aAMPs, aISIs, eps)
        %     set(gca, 'Yscale', 'log')
        %
        % end
        %
        % %% all fibers, average over pCa
        % iF = 2;
        % close all
        %
        % for i = iFs
        %     eps = squeeze(mean(RMSD(:,:,:,i),1,'omitnan'));
        %
        %     figure(1)
        %     nexttile
        %     surf(aAMPs, aISIs, eps)
        %     set(gca, 'Yscale', 'log')
        %
        %     xlabel('AMP')
        %     ylabel('ISI')
        % end
        
        %% average over identified force levels
        load([filename, '_SRS.mat'],'F0')
        
        % average
        for k = 1:length(iFs)
            for i = 1:length(th)-1
                id = F0(:,1,7,k) > th(i) & F0(:,1,7,k) <= th(i+1);
                
                RMSDs(i,:,:,k, p, ii) = mean(RMSD(id,:,:,k,p),1,'omitnan');
                F0s(i,:,:,k, ii) = mean(F0(id,:,:,k), 1,'omitnan');
            end
        end
    end
    
    %% correct for individual offsets
    % correct for individual offset
    for ii = 1:size(mcodes,1)
        for k = 1:length(iFs)
            RMSDc(:,:,:,k,p,ii) = RMSDs(:,:,:,k,p,ii) - mean(RMSDs(:,:,:,k,p,ii),'all','omitnan') + mean(RMSDs(:,:,:,:,p,ii),'all','omitnan');
        end
    end
    
    %%
end


%% make figure
ps = [4 4 4;
    5 5 5;
    7 7 7];

close all
clc

yrange = [0 .15];

color = get(gca,'colororder');
pcolors = flip(parula(7));
color = [color(2,:); pcolors(4:end-1,:);pcolors(4:end-1,:)];

figure(1)

ms = [7 6 6 5];
sym = 'shvo';

for p = 1:3
    % figure;
    for kk = 1:size(mcodes,1)
%         color = get(gca,'colororder');
        
        
        % effect of activation
        AMPid = 7;
        ISIid = 3;
        pCaid = 1:4;
        
        subplot(3,3,1 + (p-1)*3);
        
        if kk == 1
            plot(ones(1,2) * mean(F0s(pCaid(3),ISIid,AMPid,:, kk), 4, 'omitnan'), yrange, 'k--'); hold on
            
            for ii = pCaid
                plot(ones(1,2) * mean(F0s(ii,ISIid,AMPid,:, kk), 4, 'omitnan'), yrange,'k:')
            end
        end
        
        pCaid = 2:3;
        errorbar(mean(F0s(pCaid,ISIid,AMPid,:, kk), 4, 'omitnan'), mean(RMSDc(pCaid,ISIid,AMPid,:,ps(p,1), kk),4, 'omitnan'), ...
            std(RMSDc(pCaid,ISIid,AMPid,:,ps(p,1),kk),1,4, 'omitnan'), std(RMSDc(pCaid,ISIid,AMPid,:,ps(p,1),kk),1,4, 'omitnan'),...
            std(F0s(pCaid,ISIid,AMPid,:,kk),1,4, 'omitnan'), std(F0s(pCaid,ISIid,AMPid,:,kk),1,4, 'omitnan'),sym(kk), 'color', color(kk,:),'markerfacecolor', color(kk,:), 'markersize', ms(kk)); hold on
        
        pCaid = [1 4];
        errorbar(mean(F0s(pCaid,ISIid,AMPid,:, kk), 4, 'omitnan'), mean(RMSDc(pCaid,ISIid,AMPid,:,ps(p,1), kk),4, 'omitnan'), ...
            std(RMSDc(pCaid,ISIid,AMPid,:,ps(p,1),kk),1,4, 'omitnan'), std(RMSDc(pCaid,ISIid,AMPid,:,ps(p,1),kk),1,4, 'omitnan'),...
            std(F0s(pCaid,ISIid,AMPid,:,kk),1,4, 'omitnan'), std(F0s(pCaid,ISIid,AMPid,:,kk),1,4, 'omitnan'),sym(kk), 'color', color(kk,:),'markerfacecolor', [1 1 1], 'markersize', ms(kk)); hold on
        
        % effect of amplitude
        AMPid = [1, 3, 4];
        pCaid = 3;
        
        subplot(3,3,2 + (p-1)*3);
        
        if kk == 1
            plot(ones(1,2) * 0.0383, yrange, 'k-.'); hold on
        end
        
        errorbar(AMPs(AMPid), squeeze(mean(RMSDc(pCaid,ISIid,AMPid,:,ps(p,2),kk),4, 'omitnan')), ...
            squeeze(std(RMSDc(pCaid,ISIid,AMPid,:,ps(p,2),kk),1,4, 'omitnan')),sym(kk), 'color', color(kk,:),'markerfacecolor', [1 1 1], 'markersize', ms(kk)); hold on
        
        AMPid = 7;
        errorbar(AMPs(AMPid), squeeze(mean(RMSDc(pCaid,ISIid,AMPid,:,ps(p,2),kk),4, 'omitnan')), ...
            squeeze(std(RMSDc(pCaid,ISIid,AMPid,:,ps(p,2),kk),1,4, 'omitnan')),sym(kk), 'color', color(kk,:),'markerfacecolor', color(kk,:), 'markersize', ms(kk)); hold on
        
        % effect of ISI
        ISIid = 4:7;
        AMPid = 7;
        
        subplot(3,3,3 + (p-1)*3);
        
        if kk == 1
            plot(ones(1,2) * .1, yrange, 'k-.'); hold on
        end
        
        errorbar(ISIs(ISIid), squeeze(mean(RMSDc(pCaid,ISIid,AMPid,:,ps(p,3),kk),4, 'omitnan')), ...
            squeeze(std(RMSDc(pCaid,ISIid,AMPid,:,ps(p,3),kk),1,4, 'omitnan')),sym(kk), 'color', color(kk,:),'markerfacecolor', [1 1 1], 'markersize',ms(kk)); hold on
        
        ISIid = 1:3;
        errorbar(ISIs(ISIid), squeeze(mean(RMSDc(pCaid,ISIid,AMPid,:,ps(p,3),kk),4, 'omitnan')), ...
            squeeze(std(RMSDc(pCaid,ISIid,AMPid,:,ps(p,3),kk),1,4, 'omitnan')),sym(kk), 'color', color(kk,:),'markerfacecolor', color(kk,:), 'markersize', ms(kk)); hold on
        
        set(gca, 'XScale', 'log', 'Xlim', [5e-4 2e1])
        
        
    end
end


% make nice
figure(1)

titles = {'Isometric recovery', '------------------------------------------------', '------------------------------------------------', ...
    'Test stretch', '------------------------------------------------', '------------------------------------------------', ...
    'Overall', '------------------------------------------------', '------------------------------------------------'};

for j = 1:9
    subplot(3,3,j)
    
    
    title(titles{j})
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
    else
        ylabel('RMSD (F_{0})')
    end
    
    
    if j == 7
        %         title('Overall')
        xlabel('Activation (F_{0})')
    elseif j == 8
        xlabel('Amplitude (L_0)')
    elseif j == 9
        xlabel('Recovery time (s)')
    end
end

%
figure(1)
set(gcf,'units','normalized')
h = get(gcf,'position')
set(gcf,'position', [0.3536    0.2    0.2917    0.5])
