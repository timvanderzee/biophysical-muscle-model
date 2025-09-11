function[] = visualize_data_SRS()
% 
% [username, githubfolder] = get_paths();
% cd([githubfolder, '\biophysical-muscle-model\Data'])
load('SRS_data.mat', 'SRSrel', 'F0s', 'th')

%% summary plot
AMPs = [0 12 38 121 216 288 383 682]/10000;
ISIs = [1 10 100 316 1000 3160 10000]/1000;

iid = [1 3 4 5 7];
aid = [1:4, 7];

subplot(131)
plot(squeeze(F0s(:,7,1,:)), squeeze(SRSrel(:,7,1,:)), '.', 'color', [.5 .5 .5]); hold on

errorbar(mean(F0s(:,7,1,:),4,'omitnan'), mean(SRSrel(:,7,1,:),4,'omitnan'),...
    std(SRSrel(:,7,1,:),1,4,'omitnan'),std(SRSrel(:,7,1,:),1,4,'omitnan'),...
    std(F0s(:,7,1,:),1,4,'omitnan'),std(F0s(:,7,1,:),1,4,'omitnan'),'o', 'color', [.5 .5 .5]); hold on

box off
xlabel('Isometric force (F_0)')
ylabel('Relative stiffness')

xlim([0 1.05])
ylim([0 1.5])

for i = 1:length(th)
    xline(th(i),'k:')
end

subplot(132)
plot(AMPs(aid), squeeze(SRSrel(4,aid,1,:)), '.', 'color', [.5 .5 .5]); hold on
errorbar(AMPs(aid), squeeze(mean(SRSrel(4,aid,1,:), 4, 'omitnan')), squeeze(std(SRSrel(4,aid,1,:), 1, 4, 'omitnan')), 'o', 'color', [.5 .5 .5])
% xlim([1e-4 1e2])
box off

subplot(133)
semilogx(ISIs(iid), squeeze(SRSrel(4,7,iid,:)), '.', 'color', [.5 .5 .5]); hold on
errorbar(ISIs(iid), squeeze(mean(SRSrel(4,7,iid,:), 4, 'omitnan')), squeeze(std(SRSrel(4,7,iid,:), 1, 4, 'omitnan')), 'o', 'color', [.5 .5 .5])
xlim([1e-4 1e2])
box off

for j = 1:3
    subplot(1,3,j)
    ylim([0 1.5])
    yline(1,'k--')
end

% %% 3d plot for all fibers
% if ishandle(3), close(3); end
% figure(3)
% % surf(repmat(AMPs(aid)',1,5), repmat(ISIs(iid), 5, 1), mean(squeeze(SRSrel(2,aid,iid,:)),3,'omitnan')); hold on
% plot3(repmat(AMPs(aid)',1,5), repmat(ISIs(iid), 5, 1), mean(squeeze(SRSrel(2,aid,iid,:)),3,'omitnan'),'o'); hold on
% 
% set(gca,'YScale', 'log')
% 
% 
% %% 3d plot for 3 fibers
% if ishandle(3), close(3); end
% figure(3)
% 
% titles = {'Passive', 'Low', 'Intermediate', 'High', 'Maximal'};
% for i = [1, 2, 3, 5]
%     nexttile
%     surf(repmat(AMPs(:),1,7), repmat(ISIs, 8, 1), mean(squeeze(SRSrel(i,:,:,9:11)),3,'omitnan')); hold on
%     % plot3(repmat(AMPs(aid)',1,5), repmat(ISIs(iid), 5, 1), mean(squeeze(SRSrel(2,aid,iid,:)),3,'omitnan'),'o'); hold on
% 
%     set(gca,'YScale', 'log', 'CLim', [.6 1.2])
%     zlim([.6 1.2])
%     view(45, 45)
%     title(titles{i})
% %     xlabel('Amplitude')
% %     ylabel('Recovery time')
% end

end