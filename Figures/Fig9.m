clear all; close all; clc
savefig = 0;

[username, githubfolder] = get_paths();
th = [0 .07 .7 1.5];
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
        
%         SRSrel(i,:,:,k) = mean(SRS_post(id,:,:,k),1,'omitnan') ./ mean(SRS_pre(id,:,7,k),'all', 'omitnan');
        SRSrel(i,:,:,k) = mean(SRS_post(id,:,:,k),1,'omitnan') ./ mean(SRS_post(id,1,1,k),'all', 'omitnan');
        F0s(i,:,:,k) = mean(F0(id,:,:,k), 1,'omitnan');
    end
end

%%
% summary plot
eAMPs = [0 12 38 121 216 288 383 682]/10000;
eISIs = [1 10 100 316 1000 3160 10000]/1000;
aeAMPs = repmat(eAMPs, 7, 1);
aeISIs = repmat(eISIs(:), 1, 8);


% model  
modelnames = {'Hill_regular','biophysical_no_regular', 'biophysical_full_regular', 'biophysical_full_alternative'};
titles = {'Hill model', 'XB model', 'XB coop', 'XB coop + FD'};
versions = {'parms_v4', 'parms_v2d', 'parms_v2d', 'parms_v2d'};

% modelnames =  {'biophysical_full_alternative'};

for ii = 1:length(modelnames)
    modelname = modelnames{ii};
    
    cd([githubfolder, '\biophysical-muscle-model\Model output\SRS\', versions{ii}]);
load([modelname,'_SRS.mat'],'Stest', 'Scond', 'AMPs', 'iFs', 'pCas', 'ISIs', 'F0')

amAMPs = repmat(AMPs, length(ISIs), 1);
amISIs = repmat(ISIs(:), 1, length(AMPs));

% average
% iFsm = 9:11;
iFsm = [2,3,5,6,7,8,11];
iFsd = 1:11;

SRSrel_m = nan(length(th)-1,length(ISIs),length(AMPs),iFs(end));
F0s_m = nan(length(th)-1, length(ISIs),length(AMPs),iFs(end));

miid = ismember(ISIs,  eISIs);
maid = ismember(AMPs,  eAMPs);
  
for k = iFsm
    for i = 1:length(th)-1
        id = F0(:,1,7,k) > th(i) & F0(:,1,7,k) <= th(i+1);

%         SRSrel_m(i,:,:,k) = mean(Stest(id,:,:,k),1,'omitnan') ./ mean(Scond(id,1,7,k),'all', 'omitnan');
        SRSrel_m(i,:,:,k) = mean(Stest(id,:,:,k),1,'omitnan') ./ mean(Stest(id,1,1,k),'all', 'omitnan');
        F0s_m(i,:,:,k) = mean(F0(id,:,:,k), 1,'omitnan');
    end
end

% activation
k = 2;

fig = figure(1)
subplot(1,4,ii)
surf(amAMPs, amISIs, squeeze(mean(SRSrel_m(k,:,:,iFsm), 4, 'omitnan')),'edgecolor',[1 .5 .5],'linestyle','-','Facealpha',.7,'linewidth',.5); hold on
surf(aeAMPs,aeISIs, 0*ones(size(squeeze(aeISIs))),'facecolor', [1 1 1],'edgecolor', 'none'); hold on
surf(aeAMPs(1:4,7:8),aeISIs(1:4,7:8), 0*ones(size(squeeze(aeISIs(1:4,7:8)))),'facecolor', [.9 .9 .9],'edgecolor', 'none'); hold on

plot3(aeAMPs(:,[5:6, 8]), aeISIs(:,[5:6, 8]), 0*ones(size(aeAMPs(:,[5:6, 8]))), 'color', [.8 .8 .8],'linewidth',2)
plot3(aeAMPs([2 6],:)', aeISIs([2 6],:)', 0*ones(size(aeAMPs([2 6],:)))', 'color', [.8 .8 .8],'linewidth',2)

plot3(aeAMPs(:,[2:4, 7]), aeISIs(:,[2:4, 7]), 0*ones(size(aeAMPs(:,[2:4, 7]))), 'color', [.6 .6 .6],'linewidth',2)
plot3(aeAMPs([1 3:5 7],:)', aeISIs([1 3:5 7],:)', 0*ones(size(aeAMPs([1 3:5 7],:)))', 'color', [.6 .6 .6],'linewidth',2)

plot3([aeAMPs(3,7) aeAMPs(3,7)], [aeISIs(3,7) aeISIs(3,7)], [0 squeeze(mean(SRSrel(k,3,7,iFsd), 4, 'omitnan'))], '-', 'linewidth',1.5, 'color', [.6 .6 .6])
plot3([aeAMPs(1,1) aeAMPs(1,1)], [aeISIs(1,1) aeISIs(1,1)], [0 squeeze(mean(SRSrel(k,1,1,iFsd), 4, 'omitnan'))], '-', 'linewidth',1.5, 'color', [.6 .6 .6])

set(gca,'YScale','log', 'Clim', [.6 1.1], 'xtick', 0:.02:.1, 'ytick', [1e-3 1e-1 1e1]);
zlim([0 1.1])
axis([0 .07 8e-4 1.5e1 0 1.1])

plot3(aeAMPs, aeISIs, squeeze(mean(SRSrel(k,:,:,iFsd), 4, 'omitnan')),'o', 'color', [.5 .5 .5], 'markerfacecolor', [.5 .5 .5],'markersize',3); hold on

plot3(amAMPs(1,:), amISIs(1,:), squeeze(mean(SRSrel_m(k,1,:,iFsm), 4, 'omitnan')),'r-','linewidth',1.5)
plot3(amAMPs(:,1), amISIs(:,1), squeeze(mean(SRSrel_m(k,:,1,iFsm), 4, 'omitnan')),'r-','linewidth',1.5)
plot3(amAMPs(end,:), amISIs(end,:), squeeze(mean(SRSrel_m(k,end,:,iFsm), 4, 'omitnan')),'r-','linewidth',1.5)
plot3(amAMPs(:,end), amISIs(:,end), squeeze(mean(SRSrel_m(k,:,end,iFsm), 4, 'omitnan')),'r-','linewidth',1.5)

set(gca,'YScale','log')

%% compute R2
SSE = sum((mean(SRSrel_m(k,miid,maid,iFsm), 4, 'omitnan') - mean(SRSrel(k,:,:,iFsd), 4, 'omitnan')).^2,'all', 'omitnan')
SST = sum((mean(SRSrel(k,:,:,iFsd), 'all','omitnan') - mean(SRSrel(k,:,:,iFsd), 4, 'omitnan')).^2,'all', 'omitnan')

R2 = 1 - SSE./SST

set(gca,'fontsize',6)
title(titles{ii}, 'fontsize', 8)
if ii == 1
zlabel('Relative short-range stiffness','Fontsize',8)
else
    set(gca, 'zticklabel', {}, 'ZColor', 'none')
end

view(230-180,20)

text(.07, 1e0, 1.4, ['R^2 = ',num2str(round(R2, 2), 3)], 'fontsize', 6, 'horizontalalignment', 'center')


end

%%
figure(1)

subplot(141)
h = axes(fig,'visible','off');
set(h, 'Clim', [.6 1.1])
g = colorbar(h, 'location', 'WestOutside')


set(gcf,'units','centimeters','position',[10 10 19 7])
set(g,'units','centimeters','position',[5 2.5 0.2 1.7], 'fontsize', 6)

%%

       
if savefig
cd(['C:\Users\',username,'\OneDrive\9. Short-range stiffness\figures\MAT'])
figure(1)
exportgraphics(gcf,['Fig9.png'])
end