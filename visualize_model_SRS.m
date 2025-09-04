clear all; close all; clc

load('test_model_output.mat', 'Stest', 'Scond', 'AMPs', 'iFs', 'pCas', 'ISIs', 'F0')
SRSrel(:,:,:,:,2) = Stest(:,:,:,:,2)./Scond(:,:,AMPs == .0383,:,2);
F0s(:,:,:,:,2) = F0(:,:,:,:,2);

load('test_model_output_v4.mat','Stest', 'Scond', 'F0')
SRSrel(:,:,:,:,1) = Stest(:,:,:,:,1)./Scond(:,:,AMPs == .0383,:,1);
F0s(:,:,:,:,1) = F0(:,:,:,:,1);

%% plot for each fiber
% SRSrel = Stest./Scond(:,:,AMPs == .0383,:,:);
color = parula(8);

ISIid = 2;
pCaid = 5;

syms = 'osd+x*v<>phosd+x*v<>ph';

for kk = 1:2

% effect of activation
for iF = iFs
    figure(kk)
    for uu = 1:length(pCas)
        for jj = ISIid
            for ii = 1:length(AMPs)

                subplot(131);
                plot(F0s(uu,jj,ii,iF, kk), SRSrel(uu,jj,ii,iF,kk),syms(iF), 'color', color(ii,:),'markerfacecolor', color(ii,:)); hold on
                
            end
        end
    end
end

% effect of ISI
for iF = iFs
    for uu = pCaid
        for jj = 1:length(ISIs)
            for ii = 1:length(AMPs)

                subplot(132);
                semilogx(ISIs(jj), SRSrel(uu,jj,ii,iF,kk),syms(iF), 'color', color(ii,:),'markerfacecolor', color(ii,:)); hold on
                
            end
        end
    end
end

% effect of amplitude
for iF = iFs
    for uu = pCaid
        for jj = ISIid
            for ii = 1:length(AMPs)

                subplot(133);
                plot(AMPs(ii), SRSrel(uu,jj,ii,iF,kk),syms(iF), 'color', color(ii,:),'markerfacecolor', color(ii,:)); hold on
                
            end
        end
    end
end

% make nice
titles = {'Effect of activation', 'Effect of recovery', 'Effect of amplitude'};
xlabels = {'Activation', 'Recovery time', 'Amplitude'};
ylabels = 'Relative SRS';
for iF = iFs
%     figure(iF + 10)
for j = 1:3
    subplot(1,3,j)
    box off
    title(titles{j})
    xlabel(xlabels{j})
    ylabel(ylabels)
    ylim([.5 1])
end

if isfinite(F0s(pCaid,jj,ii,iF, kk))
subplot(131)
xline(F0s(pCaid,jj,ii,iF, kk), 'k--')
end
end

subplot(132)
xline(ISIs(ISIid), 'k--')

%% plot for all fibers
figure(kk+10)

% iref = find(AMPs == .0383);
% SRSrel = Stest./Scond(:,:,iref,:,:);

color = parula(8);

AMPid = 1:length(AMPs);
ISIid = 2;
pCaid = 4;

% effect of activation

for uu = 1:length(pCas)
    for jj = ISIid
        for ii = 1:length(AMPs)

            subplot(131);
            errorbar(mean(F0s(uu,jj,ii,:, kk), 4, 'omitnan'), mean(SRSrel(uu,jj,ii,:,kk),4, 'omitnan'), ...
                std(SRSrel(uu,jj,ii,:,kk),1,4, 'omitnan'), 'o', 'color', color(ii,:),'markerfacecolor', color(ii,:)); hold on

        end
    end
end

% effect of ISI
for uu = pCaid
    for jj = 1:length(ISIs)
        for ii = 1:length(AMPs)

            subplot(132);
            errorbar(ISIs(jj), mean(SRSrel(uu,jj,ii,:,kk),4, 'omitnan'), ...
                std(SRSrel(uu,jj,ii,:,kk),1,4, 'omitnan'),'o', 'color', color(ii,:),'markerfacecolor', color(ii,:)); hold on

            set(gca, 'XScale', 'log', 'Xlim', [5e-3 2e0])

        end
    end
end

% effect of amplitude
for uu = pCaid
    for jj = ISIid
        for ii = 1:length(AMPs)


            subplot(133);
            errorbar(AMPs(ii), mean(SRSrel(uu,jj,ii,:,kk),4, 'omitnan'), ...
                std(SRSrel(uu,jj,ii,:,kk),1,4, 'omitnan'),'o', 'color', color(ii,:),'markerfacecolor', color(ii,:)); hold on

        end
    end
end

% make nice
titles = {'Effect of activation', 'Effect of recovery', 'Effect of amplitude'};
xlabels = {'Activation', 'Recovery time', 'Amplitude'};
ylabels = 'Relative SRS';

for j = 1:3
    subplot(1,3,j)
    box off
    title(titles{j})
    xlabel(xlabels{j})
    ylabel(ylabels)
    ylim([.5 1.2])
    yline(1, 'k-')
    
    for ii = 1:length(AMPs)
        yline(mean(SRSrel(pCaid,ISIid,ii,:, kk),4, 'omitnan'), ':','color', color(ii,:))
    end
end

subplot(131)
xline(mean(F0s(pCaid,jj,ii,:, kk),4, 'omitnan'), 'k:')

subplot(132)
xline(ISIs(ISIid), 'k:')

end
