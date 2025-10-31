clear all; close all; clc
savefig = 0;

[username, githubfolder] = get_paths();

%% Data
eth = [0 .07 .25 .7 1.5];
tid = [3 1 7];

Cid = tid(1);
ISIid = tid(2);
AMPid = tid(3);

% [username, githubfolder] = get_paths();
% cd([githubfolder, '\biophysical-muscle-model\Data'])
load('SRS_data_v3.mat', 'SRS_pre', 'F0', 'SRS_post')

% average some pCas
iFs = 1:11;
%  iFs = [1 2 3, 5, 6, 7, 8, 10, 11];
SRSrel = nan(length(eth)-1,7,8,length(iFs));
F0s = nan(length(eth)-1, 7,8,length(iFs));

for k = iFs
    for i = 1:length(eth)-1
        id = F0(:,1,7,k) > eth(i) & F0(:,1,7,k) <= eth(i+1);
        
        SRSrel(i,:,:,k) = mean(SRS_post(id,:,:,k),1,'omitnan') ./ mean(SRS_pre(id,:,7,k),'all', 'omitnan');
        %         SRSrel(i,:,:,k) = mean(SRS_post(id,:,:,k),1,'omitnan') ./ mean(SRS_post(id,1,1,k),'all', 'omitnan');
        F0s(i,:,:,k) = mean(F0(id,:,:,k), 1,'omitnan');
    end
end

eAMPs = [0 12 38 121 216 288 383 682]/10000;
eISIs = [1 10 100 316 1000 3160 10000]/1000;

%% Model
th = [0 .07 .15 .3 .5 .7 1.5];
sACTis = [1 2 4 6];
filenames = {'Hill_regular_SRS', 'biophysical_no_regular_SRS', 'biophysical_full_regular_SRS', 'biophysical_full_alternative_SRS'};
versions = {'parms_v4', 'parms_v2d', 'parms_v2d', 'parms_v2d'};

for jj = 1:2
    if jj == 1
versions = {'parms_v4', 'parms_v4', 'parms_v4', 'parms_v4'};
    else
versions = {'parms_v4', 'parms_v2d', 'parms_v2d', 'parms_v2d'};
    end
    
for kk = 1:length(filenames)
    
    cd(['C:\Users\u0167448\Documents\GitHub\biophysical-muscle-model\Model output\SRS\', versions{kk}])
    load(filenames{kk},'Stest', 'Scond', 'AMPs', 'pCas', 'ISIs', 'F0')
    
    % average
    SRSrel_m = nan(length(th)-1,length(ISIs),length(AMPs),iFs(end));
    F0s_m = nan(length(th)-1, length(ISIs),length(AMPs),iFs(end));
    
    for k = iFs
        for i = 1:length(th)-1
            id = F0(:,1,7,k) > th(i) & F0(:,1,7,k) <= th(i+1);
            
            SRSrel_m(i,:,:,k) = mean(Stest(id,:,:,k),1,'omitnan') ./ mean(Scond(id,1,7,k),'all', 'omitnan');
            F0s_m(i,:,:,k) = mean(F0(id,:,:,k), 1,'omitnan');
        end
    end
    
    for iii = 1:3

        % variable of interest
        y = SRSrel;
        
        if iii == 1 % HD
            sAMPs = 0.0383;
            sISIs = [0.001, 0.1, 0.3160];    
        else % no HD or all
            sAMPs = [0, 0.0012, 0.0038, 0.0121, 0.0383];
            sISIs = [0.001, 0.1, 0.3160, 1, 10];
        end
                  
        miid = ismember(ISIs,  sISIs);
        maid = ismember(AMPs,  sAMPs);
        
        eiid = ismember(eISIs, sISIs);
        eaid = ismember(eAMPs, sAMPs);

        if iii == 2 % no HD
            % exclude history dependent trials
            y(:,eISIs < 1, eAMPs > .03,:) = nan;
        end
            
        % over all conditions in the dataset
        mSRSrel = mean(y, 4, 'omitnan'); % average over fibers
        N(kk,iii) = sum(isfinite(mSRSrel(:,eiid,eaid)),'all','omitnan');
        mSRSrel_m = mean(SRSrel_m, 4, 'omitnan');
        
        % average over conditions
        mmSRSrel = mean(mSRSrel(:,eiid,eaid), 'all', 'omitnan');
        
        oSST(1,iii) = sum((mSRSrel(:,eiid,eaid) - mmSRSrel).^2,'all','omitnan');
        oSSE(kk,iii) = sum((mSRSrel(:,eiid,eaid) - mSRSrel_m(sACTis,miid,maid)).^2,'all','omitnan');
        
    end
    
    % for each trial individually
    oSSEs(:,:,:,kk) = (mSRSrel(:,eiid,eaid) - mSRSrel_m(sACTis,miid,maid)).^2;
    
    
    atitles = {'Passive', 'Low activation', 'Medium activation', 'Max. activation'};
    
    for iii = 1:4
        figure(1)
        subplot(4,4,iii + (kk-1)*4)
        surf(repmat(sAMPs, length(sISIs), 1)*100, repmat(sISIs(:), 1, length(sISIs)), squeeze((mSRSrel(iii,eiid,eaid) - mSRSrel_m(sACTis(iii),miid,maid)).^2))
        set(gca,'YScale','log', 'Clim', [0 .1])
        zlim([0 .2])
        
        yt = get(gca,'Ytick');
        set(gca,'Yticklabel',yt);
        
        if kk == 1
            title(atitles{iii})
        end
        
        if iii > 1
            set(gca,'Zticklabel', [])
        end
        
        if kk < 4
            set(gca,'Xticklabel', [], 'Yticklabel', [])
        end
    end
end


%% overal R2 and AIC
id = 3;
clc
R2 = 1 - oSSE(:,id) ./ oSST(id);

% n = numel(mSRSrel_m);
n = N(1,id);

k = [5 8 9 11]';

% AIC = 2*k + n.*log(oSSE(:,id)./oSST(id))
AIC = 2*k + n.*log(oSSE(:,id)./n)

dAIC(:,jj) = AIC - AIC(1,:)

moSSE(:,:,jj) = oSSE ./ N

end

return
%% visualize

close all
figure(2)
   

color = get(gca,'colororder');
pcolors = flip(parula(7));
color = [color(2,:); pcolors(4:end-1,:);pcolors(4:end-1,:)];

bw = .25;
titles =  {'History-dependent trials', 'History-independent trials', 'All trials'};

js = [3 2 1];


for j = 1:3
    subplot(1,3,j)
    
    for i = 1:4
        bar(i, moSSE(i,js(j),2)', bw, 'FaceColor', color(i,:), 'Edgecolor', color(i,:)); hold on
    end
    
end

for j = 1:3
    subplot(1,3,j)

    for i = 1:4
        if i > 1
            fcolor = [1 1 1];
        else
            fcolor = color(i,:);
        end
        
        bar(i+bw, moSSE(i,js(j),1)', bw, 'Edgecolor', color(i,:), 'FaceColor', fcolor); hold on
    end
    
     set(gca,'fontsize',6)
     
     title(titles{js(j)}, 'fontsize', 8)
     set(gca, 'xticklabel', {})
     
     if j == 1
         ylabel('Short-range-stiffness error', 'fontsize', 8)
     else
         set(gca, 'yticklabel', {}, 'yColor', 'none')
     end
     
%     xlabel('Model', 'fontsize', 8)
    box off
    ylim([0 .12])
end

subplot(131)
% legend('Hill model','XB model', 'XB coop', 'XB coop + FD', 'location', 'best', 'fontsize', 8)
% legend boxoff

%
models = {'Hill', 'XB', 'XB coop', 'XB coop + FD'};

figure(2)
for i = 1:4
    h1(i) = rectangle('Position', [2 .12-i/60 1 .01], 'facecolor', color(i,:), 'edgecolor', color(i,:))
    
    text(3.2,  .12-i/60 + .006, models{i}, 'fontsize', 6);

    
    if i > 1
        h2(i) = rectangle('Position', [2 .12-i/60 1 .005], 'facecolor', [1 1 1], 'edgecolor', color(i,:))
   
        text(1.9,  .12-i/60 + .009, 'Original', 'fontsize', 6, 'horizontalalignment', 'right');
        text(1.9,  .12-i/60 + .003, 'Approximated', 'fontsize', 6, 'horizontalalignment', 'right');
    end
end

figure(2)
set(gcf,'units','centimeters','position',[10 10 19 5])


%%
if savefig
cd(['C:\Users\',username,'\OneDrive\9. Short-range stiffness\figures\MAT'])
       
figure(2)
exportgraphics(gcf,['Fig10.png'])
end

