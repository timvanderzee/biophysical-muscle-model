clear all; close all; clc
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
AMPs = [0 12 38 121 216 288 383 682]/10000;
ISIs = [1 10 100 316 1000 3160 10000]/1000;
aAMPs = repmat(AMPs, 7, 1);
aISIs = repmat(ISIs(:), 1, 8);

% model  
modelnames = {'Hill_regular','biophysical_no_regular', 'biophysical_full_regular', 'biophysical_full_alternative'};
% modelnames =  {'biophysical_full_alternative'};

for ii = 1:length(modelnames)
    modelname = modelnames{ii};
    
load([modelname,'_SRS.mat'],'Stest', 'Scond', 'AMPs', 'iFs', 'pCas', 'ISIs', 'F0')

% average
% iFsm = 9:11;
iFsm = [2,3,5,6,7,8,11];
iFsd = 1:11;

SRSrel_m = nan(length(th)-1,7,8,iFs(end));
F0s_m = nan(length(th)-1, 7,8,iFs(end));

for k = iFsm
    for i = 1:length(th)-1
        id = F0(:,1,7,k) > th(i) & F0(:,1,7,k) <= th(i+1);

%         SRSrel_m(i,:,:,k) = mean(Stest(id,:,:,k),1,'omitnan') ./ mean(Scond(id,1,7,k),'all', 'omitnan');
        SRSrel_m(i,:,:,k) = mean(Stest(id,:,:,k),1,'omitnan') ./ mean(Stest(id,1,1,k),'all', 'omitnan');
        F0s_m(i,:,:,k) = mean(F0(id,:,:,k), 1,'omitnan');
    end
end

k = 2;
figure(1)
nexttile
s = surf(aISIs, aAMPs, squeeze(mean(SRSrel_m(k,:,:,iFsm), 4, 'omitnan')),'edgecolor',[1 0 0],'linestyle',':','Facealpha',.7,'linewidth',.5); hold on
set(gca,'XScale','log', 'Clim', [.7 1.1])
zlim([.5 1.1])

plot3(aISIs, aAMPs, squeeze(mean(SRSrel(k,:,:,iFsd), 4, 'omitnan')),'o', 'color', [.5 .5 .5], 'markerfacecolor', [.5 .5 .5],'markersize',5); hold on
plot3(aISIs(1,:), aAMPs(1,:), squeeze(mean(SRSrel_m(k,1,:,iFsm), 4, 'omitnan')),'r-','linewidth',1)
plot3(aISIs(:,1), aAMPs(:,1), squeeze(mean(SRSrel_m(k,:,1,iFsm), 4, 'omitnan')),'r-','linewidth',1)
plot3(aISIs(end,:), aAMPs(end,:), squeeze(mean(SRSrel_m(k,end,:,iFsm), 4, 'omitnan')),'r-','linewidth',1)
plot3(aISIs(:,end), aAMPs(:,end), squeeze(mean(SRSrel_m(k,:,end,iFsm), 4, 'omitnan')),'r-','linewidth',1)

set(gca,'XScale','log')

%% compute R2
SSE = sum((mean(SRSrel_m(k,:,:,iFsm), 4, 'omitnan') - mean(SRSrel(k,:,:,iFsd), 4, 'omitnan')).^2,'all', 'omitnan')
SST = sum((mean(SRSrel(k,:,:,iFsd), 'all','omitnan') - mean(SRSrel(k,:,:,iFsd), 4, 'omitnan')).^2,'all', 'omitnan')

R2 = 1 - SSE./SST

set(gca,'fontsize',6)
xlabel('Recovery time (s)','Fontsize',8)
ylabel('Amplitude (L_0)','Fontsize',8)
zlabel('Relative short-range stiffness','Fontsize',8)
view(230,20)

end

%%
figure(1)

set(gcf,'units','normalized','position', [.2 .2 .6 .2])

% cd(['C:\Users\',username,'\OneDrive\9. Short-range stiffness\figures\MAT'])
       
% figure(1)
% exportgraphics(gcf,['Fig9B.png'])