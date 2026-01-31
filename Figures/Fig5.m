function [] = Fig5(githubfolder)


mcodes = [2 2 1; 2 1 1; 1 1 3; 1 1 2; 1 1 1; 1 2 1];
color = get_colors();

modelnames = {'Hill (no SE)', 'Hill (with SE)', '2-state XB', '2-state coop', '3-state coop', '4-state coop'};

savefig = 1;
iFs = [2,3,5,6,7,8,11];
discretized_model = 1;

if discretized_model
    versions = {'parms', 'parms_v4','parms_v4d', 'parms_v4d', 'parms_v4d', 'parms_v4d'};
else
    versions = {'parms', 'parms_v4', 'parms_v4','parms','parms_v6', 'parms_v6'};
end

%% load RMSD
% thresholds
th = [0 .07 .25 .7 1.5];

RMSDs = nan(7,7,8,11, 7, size(mcodes,1));
F0s = nan(7, 7,8,11, size(mcodes,1));
RMSDc = nan(7,7,8,11, 7, size(mcodes,1));

for p = 1:7 % phases
    for ii = 1:size(mcodes,1)
        
        mcode = mcodes(ii,:);
        
        % models
        models = {'biophysical','Hill', 'PE'};

        % model variations
        mvars = {'regular','alternative'};

        % type of cooperativity
        cvars = {'full','thin','no'};

        % choose the model version
        model = models{mcode(1)};
        mvar = mvars{mcode(2)};
        cvar = cvars{mcode(3)};
        
        if strcmp(model,'biophysical')
            filename = [model,'_',cvar,'_',mvar];
        else
           filename = [model,'_',mvar];
        end

        cd([githubfolder, '\biophysical-muscle-model\Model output\RMSD\', versions{ii}])
        load([filename, '_RMSD.mat'], 'RMSD', 'AMPs', 'ISIs', 'pCas')
        
        sAMPs = AMPs;
        sISIs = ISIs;
        
        % remove some trials
        % second fiber pCa 4.5, all small AMPs or large ISIs
        RMSD(1,:,sAMPs < .03,2,:) = nan;
        RMSD(1,sISIs > 1,:,2,:) = nan;
        
        % average over identified force levels
        sISI = .001;
        sAMP = .0383;
        
        cd([githubfolder, '\biophysical-muscle-model\Model output\SRS\Miscellaneous\', versions{ii}])
        load([filename, '_SRS.mat'],'F0', 'AMPs', 'ISIs')
        
        % average
        for k = iFs
            for i = 1:size(RMSD,1)
                %                 id = F0(:,ISIs == sISI, AMPs == sAMP,k) > th(i) & F0(:,ISIs == sISI, AMPs == sAMP,k) <= th(i+1);
                
                F0s(i,:,:,k, ii) = mean(F0(i,ismember(ISIs, sISIs),ismember(AMPs, sAMPs),k), 1,'omitnan');
                
                %                 RMSDs(i,:,:,k, p, ii) = mean(RMSD(i,:,:,k,p),1,'omitnan');
                RMSDs(i,:,:,k, p, ii) = RMSD(i,:,:,k,p);
            end
        end
    end
    
    % correct for individual offsets
    for ii = 1:size(mcodes,1)
        for k = iFs
            RMSDc(:,:,:,k,p,ii) = RMSDs(:,:,:,k,p,ii) - mean(RMSDs(:,:,:,k,p,ii),'all','omitnan') + mean(RMSDs(:,:,:,:,p,ii),'all','omitnan');
        end
    end
end


%% make figure
figure(5)

ms = [4 4 6 4 5 5];
sym = {'<','>','h','v','o','s'};

ps = [7 7 7;
    2 2 2;
    5 5 5];

%% vertical lines
yrange = [-.02 .17];

% kk = 1;
% for p = 1:3
%     AMPid = 7;
%     ISIid = 1;
%     pCaid = 1:4;
%
%     subplot(3,3,1 + (p-1)*3);
%     plot(ones(1,2) * mean(F0s(pCaid(3),ISIid,AMPid,:, kk), 4, 'omitnan'), yrange, 'k--'); hold on
%
%     subplot(3,3,2 + (p-1)*3);
%     plot(ones(1,2) * 0.0383, yrange, 'k--'); hold on
%
%     subplot(3,3,3 + (p-1)*3);
%     plot(ones(1,2) * .001, yrange, 'k--'); hold on
% end

%%
% RMSDc = RMSDs;

% close all
figure(5)

y1 = .4;

dx = nan(3, 10);

for p = 1:3
    
    % effect of activation
    AMPid = 7;
    ISIid = 1;
    
    subplot(3,4,2 + (p-1)*4);
    dydx = -y1 / 5;
    bf = y1 + dydx * (1:7);
    
    for kk = 1:6
        h = bar(kk, flip(squeeze(mean(RMSDc(:,ISIid,AMPid,:,ps(p,1), kk),4, 'omitnan')))'); hold on
        
        if kk == 1
            %            dx = nan(1, length(h));
            for i = 1:length(h)
                dx(1,i) = h(i).XOffset;
            end
        end
        
        for i = 1:length(h) % loop over bars
            errorbar(kk + dx(1,i), squeeze(mean(RMSDc(end-i+1,ISIid,AMPid,:,ps(p,1), kk),4, 'omitnan')),0,squeeze(std(RMSDc(end-i+1,ISIid,AMPid,:,ps(p,1), kk),1,4, 'omitnan')), 'marker', 'none', 'CapSize',1, 'color', [.5 .5 .5])
            
            if i == 5
                set(h(i), 'Edgecolor', brighten(color(kk,:), bf(i)), 'Facecolor', [1 1 1])
            else
                set(h(i), 'Facecolor', brighten(color(kk,:), bf(i)), 'Edgecolor', 'none', 'linewidth', 2)
            end
        end
    end
    
    % effect of amplitude
    AMPid = [1, 2, 3, 4, 7];
    pCaid = 3;
    
    subplot(3,4,3 + (p-1)*4);
    dydx = -y1 / 5;
    bf = y1 + dydx * (1:length(AMPid));
    
    
    for kk = 1:6
        h = bar(kk, squeeze(mean(RMSDc(pCaid,ISIid,AMPid,:,ps(p,2),kk),4, 'omitnan'))); hold on
        
        if kk == 1
            %             dx = nan(1, length(h));
            for i = 1:length(h)
                dx(2,i) = h(i).XOffset;
            end
        end
        
        for i = 1:length(h) % loop over bars
            errorbar(kk + dx(2,i), squeeze(mean(RMSDc(pCaid,ISIid,AMPid(i),:,ps(p,1), kk),4, 'omitnan')),0,squeeze(std(RMSDc(pCaid,ISIid,AMPid(i),:,ps(p,1), kk),1,4, 'omitnan')), 'marker', 'none', 'CapSize',1, 'color', [.5 .5 .5])
            
            if i == 5
                set(h(i), 'Edgecolor', brighten(color(kk,:), bf(i)), 'Facecolor', [1 1 1])
            else
                set(h(i), 'Facecolor', brighten(color(kk,:), bf(i)), 'Edgecolor', 'none', 'linewidth', 2)
            end
        end
    end
    
    % effect of ISI
    ISIid = [1, 3:5, 7];
    AMPid = 7;
    
    subplot(3,4,4 + (p-1)*4);
    dydx = -y1 / 2;
    bf = y1 + dydx * (1:length(ISIid));
    
    for kk = 1:6
        h = bar(kk, squeeze(mean(RMSDc(pCaid,ISIid,AMPid,:,ps(p,3),kk),4, 'omitnan'))); hold on
        
        if kk == 1
            %             dx = nan(1, length(h));
            for i = 1:length(h)
                dx(3,i) = h(i).XOffset;
            end
        end
        
        for i = 1:length(h) % loop over bars
            errorbar(kk + dx(3,i), squeeze(mean(RMSDc(pCaid,ISIid(i),AMPid,:,ps(p,1), kk),4, 'omitnan')),0,squeeze(std(RMSDc(pCaid,ISIid(i),AMPid,:,ps(p,1), kk),1,4, 'omitnan')), 'marker', 'none', 'CapSize',1, 'color', [.5 .5 .5])
            
            if i == 1
                set(h(i), 'Edgecolor', brighten(color(kk,:), bf(i)), 'Facecolor', [1 1 1])
            else
                set(h(i), 'Facecolor', brighten(color(kk,:), bf(i)), 'Edgecolor', 'none', 'linewidth', 2)
            end
        end
    end
    
    % average over all conditions
    subplot(3,4,1 + (p-1)*4);
    
    for kk = 1:6
        h = bar(kk, mean(RMSDc(:,:,:,:,ps(p,3),kk), 'all', 'omitnan'), 'facecolor', color(kk,:), 'edgecolor', 'none'); hold on
        %         errorbar(kk, mean(RMSDc(:,:,:,:,ps(p,3),kk), 'all', 'omitnan'), 0, std(RMSDc(:,:,:,:,ps(p,3),kk),1, 'all', 'omitnan'), 'marker', 'none', 'CapSize',1, 'color', [.5 .5 .5])
    end
    
end

%%

for j = 1:12
    subplot(3,4,j)
    ylim([0 .15])
    box off
    set(gca, 'fontsize', 4)
end

ABC = 'ABCD';
for p = 1:4
    subplot(3,4,p)
    text(-1, .17, ABC(p), 'fontsize', 16)
end

for p = 1:3
    subplot(3,4,1 + (p-1)*4);
    set(gca, 'Fontsize', 6)
    set(gca, 'Xtick', [])
end

set(gca, 'Xtick', 1:6, 'XTickLabelRotation', 90)
xticklabels(gca, modelnames)

figure(5)
set(gcf,'units','centimeters','position',[10 5 18 12])

%%
x = nan(3, 6, 3);

for p = 1:3
    if p == 1
        DX = dx(p,[1,4,7]);
    else
        DX = dx(p,[1,3,5]);
    end
    
    
    
    for kk = 1:6
        for i = 1:length(DX)
            x(i,kk, p) = kk + DX(i);
        end
    end
end

%%
AMPid = [1, 2, 3, 4, 7];

for j = 1:3
    if j == 1
        y = repmat(flip(pCas([1,4,7])), 1, 6);
    elseif j == 2
        y = repmat(sAMPs(AMPid([1,3,5])), 1, 6) * 100;
        
    elseif j == 3
        
        y = repmat(sISIs(ISIid([1,3,5])), 1, 6);
    end
    
    for i = 1:numel(x(:,:,j))
        g(j).xlabs{i} = num2str(y(i));
    end
end

%%
xlabels = {'Activation (pCa)', 'Amplitude (%L_0)', 'Recovery time (s)'};
figure(5)
for j = 1:3
    for p = 1:3
        X = x(:,:,p);
        subplot(3,4,j+1 + (p-1)*4);
        set(gca, 'Xtick', X(:), 'fontsize', 4)
        set(gca, 'Xticklabel', g(j).xlabs);
    end
    
    xlabel(xlabels{j}, 'fontsize', 8)
end

%% get mean RMSD
for kk = 1:6
    sigma(kk) = mean(RMSDc(:,:,:,:,ps(p,3),kk), 'all', 'omitnan');
end

%% quick stats
pval = nan(3, 5);
h = nan(3, 5);

for p = 1:3
    
    % difference between subsequent models
    dRMSDc = diff(RMSDc(:,:,:,:,ps(p,3),:), 1, 6);
    
    for i = 1:size(dRMSDc,6)
        y = dRMSDc(:,:,:,:,1,i);
        [h(p,i),pval(p,i)] = ttest(y(:), zeros(size(y(:))), 'alpha' , .05/15);
    end
end


%% indicate significant results in graph
if discretized_model
    for p = 1:3
        
        
        subplot(3,4,1 + (p-1)*4);
        
        for i = 1:5
            if h(p,i)
                plot([i i+1] + [.1 -.1], [.14 .14], 'k-', 'linewidth', 1)
                plot([i i] + .1, [mean(RMSDc(:,:,:,:,ps(p,3),i), 'all', 'omitnan')+.01 .14], 'k-', 'linewidth', 1)
                plot([i+1 i+1] -.1, [mean(RMSDc(:,:,:,:,ps(p,3),i+1), 'all', 'omitnan')+.01 .14], 'k-', 'linewidth', 1)
                
                plot(i+.5, .145, 'k*', 'markersize', 4)
                
            end
        end
    end
end


%% make nice
figure(5)

for j = 1:12
    subplot(3,4,j)
    
    
    box off
    
    %     if j < 9
    %         set(gca,'Xticklabel', [])
    %     end
    
    if j ~= 1 && j ~= 5 && j ~= 9
        set(gca,'Yticklabel', [])
        set(gca, 'Ycolor', 'none')
    else
        ylabel('Force RMSD \sigma_F (F_{0})',  'Fontsize', 8)
        %         xlim([0 1.02])
    end
    
end

%% titles
subtitles = {'Entire protocol', 'Conditioning stretch', 'Test stretch'};
titles = {'All trials', 'Effect of activation', 'Effect of amplitude', 'Effect of recovery'};


for j = 1:3
    for i = 1:4
        subplot(3,4,i + (j-1)*4)
        
        if j == 1
            title(titles{i}, 'fontsize', 8)
        else
            title('')
        end
        
        if i == 1
            title(subtitles{j}, 'fontsize', 8)
        else
            %             subtitle(' ')
        end
    end
end


end

