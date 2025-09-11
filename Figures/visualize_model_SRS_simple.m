function[] = visualize_model_SRS_simple()

% load('test_model_output.mat', 'Stest', 'Scond', 'AMPs', 'iFs', 'pCas', 'ISIs', 'F0')
% SRSrel(:,:,:,:,2) = Stest(:,:,:,:,2)./Scond(:,:,AMPs == .0383,:,2);
% F0s(:,:,:,:,2) = F0(:,:,:,:,2);
% 
% load('test_model_output_v4.mat','Stest', 'Scond', 'F0')
% SRSrel(:,:,:,:,1) = Stest(:,:,:,:,1)./Scond(:,:,AMPs == .0383,:,1);
% F0s(:,:,:,:,1) = F0(:,:,:,:,1);

load('test_model_output_v6.mat','Stest', 'Scond', 'AMPs', 'iFs', 'pCas', 'ISIs', 'F0')
SRSrel = Stest./Scond(:,:,AMPs == .0383,:,2);
F0s = F0;

%% plot for each fiber
% SRSrel = Stest./Scond(:,:,AMPs == .0383,:,:);
color = parula(8);

ISIid = 2;
pCaid = 5;

syms = 'osd+x*v<>phosd+x*v<>ph';
kk = 2;


%% plot for all fibers
% figure(1)
color = parula(8);

AMPid = 7;
ISIid = 1;
pCaid = 4;

% AMPid = 5;

% effect of activation

for uu = 1:length(pCas)
    for jj = ISIid
        for ii = AMPid

            subplot(131);
            errorbar(mean(F0s(uu,jj,ii,:, kk), 4, 'omitnan'), mean(SRSrel(uu,jj,ii,:,kk),4, 'omitnan'), ...
                std(SRSrel(uu,jj,ii,:,kk),1,4, 'omitnan'), 'o', 'color', color(ii,:),'markerfacecolor', color(ii,:)); hold on

        end
    end
end


% effect of amplitude
for uu = pCaid
    for jj = ISIid
        for ii = 1:length(AMPs)

            subplot(132);
            errorbar(AMPs(ii), mean(SRSrel(uu,jj,ii,:,kk),4, 'omitnan'), ...
                std(SRSrel(uu,jj,ii,:,kk),1,4, 'omitnan'),'o', 'color', color(ii,:),'markerfacecolor', color(ii,:)); hold on

        end
    end
end

% effect of ISI
for uu = pCaid
    for jj = 1:length(ISIs)
        for ii = AMPid

            subplot(133);
            errorbar(ISIs(jj), mean(SRSrel(uu,jj,ii,:,kk),4, 'omitnan'), ...
                std(SRSrel(uu,jj,ii,:,kk),1,4, 'omitnan'),'o', 'color', color(ii,:),'markerfacecolor', color(ii,:)); hold on

            set(gca, 'XScale', 'log', 'Xlim', [5e-3 2e0])

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
    ylim([.4 1.2])
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
