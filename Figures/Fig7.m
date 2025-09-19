clear all; close all; clc

mcodes = [1 1 1; 2 1 1];

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

    
    
%     figure;

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
for ii = 1:2
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
  figure(1)
for p = 1:3
% figure;
for kk = 1:2
color = get(gca,'colororder');

AMPid = 7;
ISIid = 3;
pCaid = 3;

% effect of activation
for uu = 1:size(RMSDc, 1)
    for jj = ISIid
        for ii = AMPid

            subplot(3,3,1 + (p-1)*3);
            errorbar(mean(F0s(uu,jj,ii,:, kk), 4, 'omitnan'), mean(RMSDc(uu,jj,ii,:,ps(p,1), kk),4, 'omitnan'), ...
                std(RMSDc(uu,jj,ii,:,ps(p,1),kk),1,4, 'omitnan'), std(RMSDc(uu,jj,ii,:,ps(p,1),kk),1,4, 'omitnan'),...
                std(F0s(uu,jj,ii,:,kk),1,4, 'omitnan'), std(F0s(uu,jj,ii,:,kk),1,4, 'omitnan'),'o', 'color', color(kk,:),'markerfacecolor', color(kk,:)); hold on

            if uu == 3
                    xline(mean(F0s(uu,jj,ii,:, kk), 4, 'omitnan'), 'k--')
                    
            end
            xline(mean(F0s(uu,jj,ii,:, kk), 4, 'omitnan'), 'k:')
                    
        end
    end
end

% effect of amplitude
for uu = pCaid
    for jj = ISIid
        for ii = [1, 3, 4, 7]

            subplot(3,3,2 + (p-1)*3);
            errorbar(AMPs(ii), mean(RMSDc(uu,jj,ii,:,ps(p,2),kk),4, 'omitnan'), ...
                std(RMSDc(uu,jj,ii,:,ps(p,2),kk),1,4, 'omitnan'),'o', 'color', color(kk,:),'markerfacecolor', color(kk,:)); hold on

             xline(0.0383, 'k--')
             xline(0.0383, 'k:')
        
        end
    end
end

% effect of ISI
for uu = pCaid
    for jj = 1:length(ISIs)
        for ii = AMPid

            subplot(3,3,3 + (p-1)*3);
            errorbar(ISIs(jj), mean(RMSDc(uu,jj,ii,:,ps(p,3),kk),4, 'omitnan'), ...
                std(RMSDc(uu,jj,ii,:,ps(p,3),kk),1,4, 'omitnan'),'o', 'color', color(kk,:),'markerfacecolor', color(kk,:)); hold on

            set(gca, 'XScale', 'log', 'Xlim', [5e-4 2e1])
            xline(0.1, 'k--')
            xline(0.1, 'k:')
    
        end
    end
end
end


% % make nice
% titles = {'Effect of activation', 'Effect of amplitude','Effect of recovery'};
% xlabels = {'Activation', 'Amplitude', 'Recovery time'};
% ylabels = 'RMSD';



% %     yline(1, 'k-')
%     
% %     for ii = 1:length(AMPs)
% %         yline(mean(SRSrel(pCaid,ISIid,ii,:, kk),4, 'omitnan'), ':','color', color(ii,:))
% %     end
% end
% 
% subplot(131)
% % xline(mean(F0s(pCaid,jj,ii,:, kk),4, 'omitnan'), 'k:')
% 
% subplot(132)
% % xline(AMPs(AMPid), 'k:')
% 
% subplot(133)
% % xline(ISIs(ISIid), 'k:')
% 
% xlim([1e-4 1e1])

end


%% make nice
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
    ylim([0 .15])
    
    
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


