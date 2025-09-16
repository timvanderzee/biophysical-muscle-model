function[] = visualize_model_SRS_simple()

% load('test_model_output.mat', 'Stest', 'Scond', 'AMPs', 'iFs', 'pCas', 'ISIs', 'F0')
% SRSrel(:,:,:,:,2) = Stest(:,:,:,:,2)./Scond(:,:,AMPs == .0383,:,2);
% F0s(:,:,:,:,2) = F0(:,:,:,:,2);
% 
% load('test_model_output_v4.mat','Stest', 'Scond', 'F0')
% SRSrel(:,:,:,:,1) = Stest(:,:,:,:,1)./Scond(:,:,AMPs == .0383,:,1);
% F0s(:,:,:,:,1) = F0(:,:,:,:,1);

filenames = {'test_model_output_v8.mat', 'Hill_regular_SRS'};

for ff = 1:length(filenames)

load(filenames{ff},'Stest', 'Scond', 'AMPs', 'iFs', 'pCas', 'ISIs', 'F0')
% SRSrel = Stest./Scond(:,:,AMPs == .0383,:);
% F0s = F0;

% average 
iFs = [1 2 3, 5, 6, 7, 8, 10, 11];
th = [0 .05 .1 .25 .7 1.5];
SRSrel = nan(length(th)-1,7,8,length(iFs));
F0s = nan(length(th)-1, 7,8,length(iFs));

for k = 1:length(iFs)
    for i = 1:length(th)-1
        id = F0(:,1,7,k) > th(i) & F0(:,1,7,k) <= th(i+1);
        
        SRSrel(i,:,:,k) = mean(Stest(id,:,:,k),1,'omitnan') ./ mean(Scond(id,1,7,k),'all', 'omitnan');
        F0s(i,:,:,k) = mean(F0(id,:,:,k), 1,'omitnan');
    end
end

%% plot for all fibers
kk = 1;
% figure(1)
% color = parula(8);
color = get(gca,'colororder');

AMPid = 7;
ISIid = 1;
pCaid = 3;

% subplot(131)
% for uu = 1:size(SRSrel, 1)
%     for jj = ISIid
%         for ii = AMPid
%             for i = 1:size(SRSrel, 4)
%             
%                 plot(F0s(:,jj,ii,i, kk), SRSrel(:,jj,ii,i,kk),'+'); hold on
%             end
%         end
%     end
% end

% keyboard
% effect of activation
for uu = 1:size(SRSrel, 1)
    for jj = ISIid
        for ii = AMPid

            subplot(131);
            errorbar(mean(F0s(uu,jj,ii,:, kk), 4, 'omitnan'), mean(SRSrel(uu,jj,ii,:,kk),4, 'omitnan'), ...
                std(SRSrel(uu,jj,ii,:,kk),1,4, 'omitnan'), std(SRSrel(uu,jj,ii,:,kk),1,4, 'omitnan'),...
                std(F0s(uu,jj,ii,:,kk),1,4, 'omitnan'), std(F0s(uu,jj,ii,:,kk),1,4, 'omitnan'),'o', 'color', color(ff,:),'markerfacecolor', color(ff,:)); hold on

        end
    end
end


% effect of amplitude
for uu = pCaid
    for jj = ISIid
        for ii = 1:length(AMPs)

            subplot(132);
            errorbar(AMPs(ii), mean(SRSrel(uu,jj,ii,:,kk),4, 'omitnan'), ...
                std(SRSrel(uu,jj,ii,:,kk),1,4, 'omitnan'),'o', 'color', color(ff,:),'markerfacecolor', color(ff,:)); hold on

        end
    end
end

% effect of ISI
for uu = pCaid
    for jj = 1:length(ISIs)
        for ii = AMPid

            subplot(133);
            errorbar(ISIs(jj), mean(SRSrel(uu,jj,ii,:,kk),4, 'omitnan'), ...
                std(SRSrel(uu,jj,ii,:,kk),1,4, 'omitnan'),'o', 'color', color(ff,:),'markerfacecolor', color(ff,:)); hold on

            set(gca, 'XScale', 'log', 'Xlim', [5e-3 2e0])

        end
    end
end
end

% make nice
titles = {'Effect of activation', 'Effect of amplitude','Effect of recovery'};
xlabels = {'Activation', 'Amplitude', 'Recovery time'};
ylabels = 'Relative SRS';

for j = 1:3
    subplot(1,3,j)
    box off
    title(titles{j})
    xlabel(xlabels{j})
    ylabel(ylabels)
    ylim([0 2])
    yline(1, 'k-')
    
%     for ii = 1:length(AMPs)
%         yline(mean(SRSrel(pCaid,ISIid,ii,:, kk),4, 'omitnan'), ':','color', color(ii,:))
%     end
end

subplot(131)
% xline(mean(F0s(pCaid,jj,ii,:, kk),4, 'omitnan'), 'k:')

subplot(132)
% xline(AMPs(AMPid), 'k:')

subplot(133)
% xline(ISIs(ISIid), 'k:')

xlim([1e-4 1e1])

end
