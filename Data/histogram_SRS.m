clear all; close all; clc

load('SRS_data.mat')

%%
SRS = SRSrel(3,:,:,:);

AMPs = [0 12 38 121 216 288 383 682]/10000;
ISIs = [1 10 100 316 1000 3160 10000]/1000;

xAMPs = repmat(AMPs, length(ISIs),1);
yISIs = repmat(ISIs(:), 1, length(AMPs),1);

figure(1)
surf(xAMPs, yISIs, squeeze(mean(SRSrel(4,:,:,:), 4, 'omitnan')))

set(gca,'YScale','log')

%%
close all
figure(2)

% for i = 3:(size(SRSrel,1))
%     nexttile
SRS = SRSrel(3:end,:,:,:);
histogram(SRS(:),'BinWidth',.05); hold on

SRS2 = SRSrel(3:end,ISIs < 1, AMPs > .03,:);
histogram(SRS2,'BinWidth',.05); hold on
xlim([.3 1.5])
box off
xlabel('Relative short-range stiffness')
ylabel('Count')
legend('All trials', '(ISI < 1 s) AND (AMP > 3% L_0)', 'location','best')
legend boxoff
% end