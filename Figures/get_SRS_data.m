clear all; close all; clc

usernames = {'timvd','u0167448'};

for i = 1:length(usernames)
    docfolder =  ['C:\Users\', usernames{i}, '\Documents'];
    
    if isfolder(docfolder)
        mainfolder = docfolder;
        username = usernames{i};
    end
end

if contains(mainfolder, 'timvd')
    githubfolder = mainfolder;
else
    githubfolder = [mainfolder, '\GitHub'];
end

addpath(genpath([githubfolder, '\muscle-thixotropy']))
addpath(genpath([mainfolder, '\casadi-3.7.1-windows64-matlab2018b']))
addpath(genpath([githubfolder, '\biophysical-muscle-model']))

%%
Kss = 1:7;
tiso = 3;

iFs = 1:11;

F0 = nan(length(Kss), 8,7,length(iFs));
SRS_pre = nan(length(Kss), 8,7,length(iFs));
SRS_post = nan(length(Kss), 8,7,length(iFs));
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

for k = 1:length(iFs)
    cd(['C:\Users\',username,'\OneDrive - KU Leuven\9. Short-range stiffness\matlab\data'])
    load([fibers{iFs(k)},'_cor_new.mat'],'data')
    
    disp([fibers{iFs(k)}])
    
    for n = 1:7
        disp(n)
        for m = 1:8
            
            Data = prep_data(data,n,m,Kss,tiso);
            
            if ~isempty(Data.t)
                ts = 0:tiso:(tiso*(length(Kss)-1));
                
                % pre-allocate
                ns = nan(length(Kss), 2);
                os = nan(length(Kss), 2);
                ds = nan(length(Kss), 2);
                
                [id0,id1,id2] = get_indices(Data.t, tiso, ts, Data.dTt, Data.dTc, Data.ISI, Data.Ca(Kss));
                
%                 id0 = repmat(1:100,7,1);
%                     figure(1)
%                     nexttile
%                     plot(Data.L, Data.F,'.'); hold on
                
                for i = 1:length(Kss)
                    
                    % filter
                    fs = 1000;
                    fc = 50;
                    Wn = fc / (.5*fs);
                    [b,a] = butter(2, Wn);
                    
                    F = nan(size(Data.F));
                    L = nan(size(Data.F));
                    
                    id = isfinite(Data.F);
                    F(id) = filtfilt(b,a, Data.F(id));
                    L(id) = filtfilt(b,a, Data.L(id));
                    
                    dp1 = polyfit(L(id1(i,:)), F(id1(i,:)), 1);
                    
                    SRS_post(i,m,n,k) = dp1(1);
%                     plot(L(id1(i,:)), F(id1(i,:)), 'r.')
                    
                    
                    % only for ISI = 1 ms
                    if n == 1
                        dp2 = polyfit(L(id2(i,:)), F(id2(i,:)), 1);
                        
                        F0(i,m,n,k) = mean(F(id0(i,:)));
                        SRS_pre(i,m,n,k) = dp2(1);
   
%                         plot(L(id0(i,:)), F(id0(i,:)), 'g.')
%                         plot(L(id2(i,:)), F(id2(i,:)), 'y.')
                    end
                end
            end
        end
    end
end

%% average some pCas
th = [0 .05 .1 .25 .7 1.5];
SRSrel = nan(length(th)-1,8,7,length(iFs));
F0s = nan(length(th)-1, 8,7,length(iFs));

for k = 1:length(iFs)
    for i = 1:length(th)-1
        id = F0(:,7,1,k) > th(i) & F0(:,7,1,k) <= th(i+1);
        
        SRSrel(i,:,:,k) = mean(SRS_post(id,:,:,k),1,'omitnan') ./ mean(SRS_pre(id,:,1,k),'all', 'omitnan');
        F0s(i,:,:,k) = mean(F0(id,:,:,k), 1,'omitnan');
    end
end


%% summary plot
AMPs = [0 12 38 121 216 288 383 682]/10000;
ISIs = [1 10 100 316 1000 3160 10000]/1000;

iid = [1 3 4 5 7];
aid = [1:4, 7];

if ishandle(2), close(2); end
figure(2)

subplot(131)
plot(squeeze(F0s(:,7,1,:)), squeeze(SRSrel(:,7,1,:)), '.', 'color', [.5 .5 .5]); hold on

errorbar(mean(F0s(:,7,1,:),4,'omitnan'), mean(SRSrel(:,7,1,:),4,'omitnan'),...
    std(SRSrel(:,7,1,:),1,4,'omitnan'),std(SRSrel(:,7,1,:),1,4,'omitnan'),...
    std(F0s(:,7,1,:),1,4,'omitnan'),std(F0s(:,7,1,:),1,4,'omitnan'),'o-'); hold on

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
errorbar(AMPs(aid), squeeze(mean(SRSrel(4,aid,1,:), 4, 'omitnan')), squeeze(std(SRSrel(4,aid,1,:), 1, 4, 'omitnan')), 'o-')
% xlim([1e-4 1e2])
box off

subplot(133)
semilogx(ISIs(iid), squeeze(SRSrel(4,7,iid,:)), '.', 'color', [.5 .5 .5]); hold on
errorbar(ISIs(iid), squeeze(mean(SRSrel(4,7,iid,:), 4, 'omitnan')), squeeze(std(SRSrel(4,7,iid,:), 1, 4, 'omitnan')), 'o-')
xlim([1e-4 1e2])
box off

for j = 1:3
    subplot(1,3,j)
    ylim([0 1.5])
    yline(1,'k--')
end

%% 3d plot for all fibers
if ishandle(3), close(3); end
figure(3)
% surf(repmat(AMPs(aid)',1,5), repmat(ISIs(iid), 5, 1), mean(squeeze(SRSrel(2,aid,iid,:)),3,'omitnan')); hold on
plot3(repmat(AMPs(aid)',1,5), repmat(ISIs(iid), 5, 1), mean(squeeze(SRSrel(2,aid,iid,:)),3,'omitnan'),'o'); hold on

set(gca,'YScale', 'log')


%% 3d plot for 3 fibers
if ishandle(3), close(3); end
figure(3)

titles = {'Passive', 'Low', 'Intermediate', 'High', 'Maximal'};
for i = [1, 2, 3, 5]
    nexttile
    surf(repmat(AMPs(:),1,7), repmat(ISIs, 8, 1), mean(squeeze(SRSrel(i,:,:,9:11)),3,'omitnan')); hold on
    % plot3(repmat(AMPs(aid)',1,5), repmat(ISIs(iid), 5, 1), mean(squeeze(SRSrel(2,aid,iid,:)),3,'omitnan'),'o'); hold on

    set(gca,'YScale', 'log', 'CLim', [.6 1.2])
    zlim([.6 1.2])
    view(45, 45)
    title(titles{i})
%     xlabel('Amplitude')
%     ylabel('Recovery time')
end