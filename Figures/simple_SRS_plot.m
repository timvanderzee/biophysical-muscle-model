clear all; close all; clc
savefig = 0;

[username, githubfolder] = get_paths();
color = lines(7);

% Model
tid = [3 1 7];

showbar = 0;
showline = 1;
% figure(1)
filenames = {'biophysical_full_regular_SRS', 'biophysical_full_regular_SRS_old'};
versions = {'parms_v4', 'parms_v6'};

iFs = [2,3,5,6,7,8,11];
yrange = [.3 1.5];

kk = 1;
 cd(['C:\Users\u0167448\Documents\GitHub\biophysical-muscle-model\Model output\SRS\', versions{kk}])
    load(filenames{kk},'Stest', 'Scond', 'AMPs', 'pCas', 'ISIs', 'F0')
   
    
% average
SRSrel_m = nan(length(pCas),length(ISIs),length(AMPs),iFs(end),2);
F0s_m = nan(length(pCas), length(ISIs),length(AMPs),iFs(end),2);

Flin = linspace(.02, .98, 100);

for kk = 1:length(filenames)
    
    cd(['C:\Users\u0167448\Documents\GitHub\biophysical-muscle-model\Model output\SRS\', versions{kk}])
    load(filenames{kk},'Stest', 'Scond', 'AMPs', 'pCas', 'ISIs', 'F0')
   
    for k = iFs
        for i = 1:size(Stest,1)
            SRSrel_m(i,:,:,k,kk) = Stest(i,:,:,k) ./ Scond(i,1,7,k);
            F0s_m(i,:,:,k,kk) = F0(i,:,:,k);
           
        end
    end 
end


load('SRS_data_v4.mat', 'SRS_pre', 'F0', 'SRS_post')

% average some pCas
iFs = [2,3,5,6,7,8,11];

SRSrel = nan(length(pCas),7,8,iFs(end));
F0s = nan(length(pCas), 7,8,iFs(end));

for k = iFs
    for i = 1:size(Stest,1)
 
        SRSrel(i,:,:,k) = mean(SRS_post(i,:,:,k),1,'omitnan') ./ mean(SRS_pre(i,:,7,k),'all', 'omitnan');
        F0s(i,:,:,k) = mean(F0(i,:,:,k), 1,'omitnan');
    end
end

%% lump data
iFs = [2,3,5,6,7,8,11];
eth = [0 .05 .3 .7 1.2];

SRSrelL = nan(length(eth)-1,7,8,iFs(end));
F0sL = nan(length(eth)-1, 7,8,iFs(end));

for k = iFs
    for i = 1:length(eth)-1
        id = F0(:,1,7,k) > eth(i) & F0(:,1,7,k) <= eth(i+1);
        
        SRSrelL(i,:,:,k) = mean(SRS_post(id,:,:,k),1,'omitnan') ./ mean(SRS_pre(id,:,7,k),'all', 'omitnan');
        F0sL(i,:,:,k) = mean(F0(id,:,:,k), 1,'omitnan');
    end
end

%% interpolate
SRSint = nan(length(Flin), size(SRSrel_m,2), size(SRSrel_m,3), size(SRSrel_m,4), size(SRSrel_m,5));

for i = 1:size(SRSrel_m,2)
    for j = 1:size(SRSrel_m,3)
        for k = 1:size(SRSrel_m,4)
            for kk = 1:size(SRSrel_m,5)
            if sum(isfinite(F0s_m(:,i,j,k))) > 1
                SRSint(:,i,j,k,kk) = interp1(F0s_m(:,i,j,k,kk), SRSrel_m(:,i,j,k,kk), Flin);
            end
            end
        end
    end
end

%%
close all
figure(1)

n = 1;
m = 7;

markers = 'so'; 

% color = parula(iFs(end));
color = lines(7);


for k = iFs
    
%     nexttile
    
    errorbar(mean(F0sL(:,n,m,:), 4, 'omitnan'), mean(SRSrelL(:,n,m,:), 4, 'omitnan'), std(SRSrelL(:,n,m,:), 1, 4, 'omitnan'), std(SRSrelL(:,n,m,:), 1, 4, 'omitnan'),std(F0sL(:,n,m,:), 1, 4, 'omitnan'),std(F0sL(:,n,m,:), 1, 4, 'omitnan'), 'o'); hold on

%     plot(F0s(:,n,m,k), SRSrel(:,n,m,k),'o', 'color', [.5 .5 .5], 'markerfacecolor', [1 1 1]); hold on
    
    for kk = 1:2
%         plot(F0s_m(:,n,m,k,kk), SRSrel_m(:,n,m,k,kk), markers(kk), 'color', color(kk,:), 'markerfacecolor', color(kk,:)); hold on
%         plot(Flin, SRSint(:,n,m,k,kk),'--', 'color', color(kk,:))
    end
 
    
    
    yline(1, 'k--')
%     ylim([0 1.5])
    box off
    ylabel('SRS (-)')
    xlabel('Activation (-)')
    title(['Fiber ', num2str(k)])
end


plot(Flin, mean(SRSint(:,n,m,:,1),4,'omitnan'), '-', 'linewidth', 2)
plot(Flin, mean(SRSint(:,n,m,:,2),4,'omitnan'), '-', 'linewidth', 2)

%% R^2
SRSi = interp1(Flin, mean(SRSint(:,n,m,:,2),4,'omitnan'), mean(F0sL(:,n,m,:), 4, 'omitnan'), [], 'extrap')


SStot = sum((mean(SRSrelL(:,n,m,:), 4, 'omitnan') - mean(SRSrelL(:,n,m,:), 'all', 'omitnan')).^2)
SSrest = sum((mean(SRSrelL(:,n,m,:), 4, 'omitnan') - SRSi).^2)

R2 = 1 - SSrest/SStot
