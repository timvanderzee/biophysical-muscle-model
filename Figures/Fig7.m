function[] = Fig7(datafolder, modelfolder, githubfolder)

% [username, githubfolder] = get_paths();
th = [0 .3 .7 1.5];
tid = [3 1 7];

%% Data
Cid = tid(1);
ISIid = tid(2);
AMPid = tid(3);

cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Data'))
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
modelnames = {'Hill_alternative', 'Hill_regular', 'biophysical_full_regular'};
titles = {'Hill-type model (no SE)', 'Hill-type model (with SE)', '3-state XB coop model'};

discretized_model = 0;
if discretized_model
    versions = {'parms', 'parms_v3','parms_v4d'};
else
    versions = {'parms', 'parms_v3', 'parms_v6'};
end

acolors = lines(6);
acolors = [acolors(end,:); acolors];

id = [1, 2, 5];
color = acolors(id,:);

for ii = 1:length(modelnames)
    modelname = modelnames{ii};
    
    if contains(modelname, 'biophysical') % biophysical model
    if discretized_model 
        subfolder =  'Discretized';
    else
        subfolder = 'Approximated';
    end
    
    else % Hill-type model
        subfolder = 'Hill';
    end
    
    cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Model output', 'SRS', subfolder));
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
    
    fig = figure(7)
    subplot(3,3,ii)
    surf(amAMPs*100, amISIs, squeeze(mean(SRSrel_m(k,:,:,iFsm), 4, 'omitnan')),'edgecolor',[.8 .8 .8],'linestyle','-','Facealpha',.7,'linewidth',.5); hold on
    surf(aeAMPs*100,aeISIs, 0*ones(size(squeeze(aeISIs))),'facecolor', [1 1 1],'edgecolor', 'none'); hold on
    surf(aeAMPs(1:4,7:8)*100,aeISIs(1:4,7:8), 0*ones(size(squeeze(aeISIs(1:4,7:8)))),'facecolor', [.9 .9 .9],'edgecolor', 'none'); hold on
    
    colors = color(ii,:) .* repmat(linspace(0,1,100)', 1, 3);
    set(gca, 'colormap', colors)
    
    plot3(aeAMPs(:,[5:6, 8])*100, aeISIs(:,[5:6, 8]), 0*ones(size(aeAMPs(:,[5:6, 8]))), 'color', [.8 .8 .8],'linewidth',2)
    plot3(aeAMPs([2 6],:)'*100, aeISIs([2 6],:)', 0*ones(size(aeAMPs([2 6],:)))', 'color', [.8 .8 .8],'linewidth',2)
    
    plot3(aeAMPs(:,[2:4, 7])*100, aeISIs(:,[2:4, 7]), 0*ones(size(aeAMPs(:,[2:4, 7]))), 'color', [.6 .6 .6],'linewidth',2)
    plot3(aeAMPs([1 3:5 7],:)'*100, aeISIs([1 3:5 7],:)', 0*ones(size(aeAMPs([1 3:5 7],:)))', 'color', [.6 .6 .6],'linewidth',2)
    
    plot3([aeAMPs(3,7) aeAMPs(3,7)]*100, [aeISIs(3,7) aeISIs(3,7)], [0 squeeze(mean(SRSrel(k,3,7,iFsd), 4, 'omitnan'))], '-', 'linewidth',1.5, 'color', [.6 .6 .6])
    plot3([aeAMPs(1,1) aeAMPs(1,1)]*100, [aeISIs(1,1) aeISIs(1,1)], [0 squeeze(mean(SRSrel(k,1,1,iFsd), 4, 'omitnan'))], '-', 'linewidth',1.5, 'color', [.6 .6 .6])
    
    set(gca,'YScale','log', 'Clim', [.6 1.1], 'xtick', (0:.02:.1)*100, 'ytick', [1e-3 1e-1 1e1]);
    zlim([0 1.1])
    axis([0 .07*100 8e-4 1.5e1 0 1.1])
    
    plot3(aeAMPs*100, aeISIs, squeeze(mean(SRSrel(k,:,:,iFsd), 4, 'omitnan')),'o', 'color', [.5 .5 .5], 'markerfacecolor', [.5 .5 .5],'markersize',3); hold on
    
    plot3(amAMPs(1,:)*100, amISIs(1,:), squeeze(mean(SRSrel_m(k,1,:,iFsm), 4, 'omitnan')),'-','linewidth',1.5, 'color', color(ii,:))
    plot3(amAMPs(:,1)*100, amISIs(:,1), squeeze(mean(SRSrel_m(k,:,1,iFsm), 4, 'omitnan')),'-','linewidth',1.5, 'color', color(ii,:))
    plot3(amAMPs(end,:)*100, amISIs(end,:), squeeze(mean(SRSrel_m(k,end,:,iFsm), 4, 'omitnan')),'-','linewidth',1.5, 'color', color(ii,:))
    plot3(amAMPs(:,end)*100, amISIs(:,end), squeeze(mean(SRSrel_m(k,:,end,iFsm), 4, 'omitnan')),'-','linewidth',1.5, 'color', color(ii,:))
    
    set(gca,'YScale','log')
    
    %% compute R2
    SSE = sum((mean(SRSrel_m(k,miid,maid,iFsm), 4, 'omitnan') - mean(SRSrel(k,:,:,iFsd), 4, 'omitnan')).^2,'all', 'omitnan')
    SST = sum((mean(SRSrel(k,:,:,iFsd), 'all','omitnan') - mean(SRSrel(k,:,:,iFsd), 4, 'omitnan')).^2,'all', 'omitnan')
    
    R2 = 1 - SSE./SST
    
    set(gca,'fontsize',6)
    subtitle(titles{ii}, 'fontsize', 8, 'color', color(ii,:))
    if ii == 1
        zlabel('Relative short-range stiffness','Fontsize',8)
    else
        set(gca, 'zticklabel', {}, 'ZColor', 'none')
    end
    
    view(230-180,20)
    
    %     text(.07, 1e0, 1.4, ['R^2 = ',num2str(round(R2, 2), 3)], 'fontsize', 6, 'horizontalalignment', 'center')
    
    
end

%%
figure(7)

subplot(331)
h = axes(fig,'visible','off');
set(h, 'Clim', [.6 1.1])
% g = colorbar(h, 'location', 'WestOutside')

% set(g,'units','centimeters','position',[6 2.5+11 0.2 1.7], 'fontsize', 6)

%% plot force traces
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

iFsm = [2,3,5,6,7,8,11];
AMP = .0383;

tlin = -2:.001:1;
Fexps = nan(length(tlin), max(iFsm),3, 2);

pCas = [4.5 6.1 6.2 6.3 6.4 6.6 9];
pCai = nan(11,4);

tiso = 2;
dTt = .0383/.4545; % test stretch (= constant);
% t0 = 2*dTt + tiso + ISI;

% eth = [0 .05 .3 .7 1.2];
thmin = [.7 .3 .05];
thmax = [2 .7 .3];

ISIs = [.001 1];



for k = 1:2
    ISI = ISIs(k);
    
    for iF = iFsm
        years = {'2017', '2018'};
        for m = 1:length(years)
            if contains(fibers{iF}, years{m})
                fullfolder = [datafolder, '\', years{m}];
            end
        end
        
        cd(fullfolder)
        load([fibers{iF},'.mat'],'data');
        
        %     figure(k+10)
        %     nexttile
        
        for j = 1:length(thmin)
            
            for i = 1:length(pCas)
                pCa = pCas(i);
                cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Figures'))
                [texp, Lexp, Fexp, Tsrel] = get_data(data, ISI, AMP, pCa);
                
                Fiso = mean(Fexp(texp<-.5), 'omitnan');
                
                if Fiso < thmax(j) && Fiso > thmin(j)
                    
                    Fexps(:,iF,j,k) = interp1(texp(isfinite(Fexp)), Fexp(isfinite(Fexp)), tlin);
                    pCai(iF,j) = pCa;
                    
                    %                 plot(texp, Fexp,'-', tlin, Fexps(:,iF,j),'.'); hold on
                    %                 xlim([-1 1])
                    break
                end
            end
        end
    end
    
    %%
    figure(7)
    
    % pCai(:,end) = 9;
    % t0 = 0;
    t0 = 2*dTt + tiso + ISI;
    
    % ISI = 1;
    oFis = nan(max(iFsm), 4000);
    
    % close all
    % figure(7)
    
    if k == 1
        mcodes = [2 2 1; 2 1 1; 1 1 1];
        id = [1, 2, 5];
        
    else
        mcodes = [1 1 1];
        id = 5;
    end
    
    color = acolors(id,:);
    
    for j = 1:size(mcodes,1)
        if k == 1
            subplot(3,3,j+3)
            xlim([-.2 .15])
            box off
            
        else
            subplot(3,3,[7 9])
            xlim([-1.2 .15])
            box off
        end
        
        hold(gca, 'on')
        
        set(gca, 'Fontsize', 6)
        xlabel(gca,'Time (s)', 'Fontsize', 8)
        
        if j == 1
            ylabel(gca, 'Force (-)', 'Fontsize', 8)
            
        end
        
        for jj = 1:3
            %     subplot(1,3,jj)
            % plot(tlin, Fexps, 'k-', 'linewidth', 1); hold on
            plot(tlin, mean(Fexps(:,:,jj,k), 2, 'omitnan'), 'k-', 'linewidth', 2); hold on
            
            
            % 2. Calculate Upper and Lower Bounds
            sd = std(Fexps(:,:,jj,k), 1,2, 'omitnan');
            upper_bound = mean(Fexps(:,:,jj,k), 2, 'omitnan') + sd;
            lower_bound = mean(Fexps(:,:,jj,k), 2, 'omitnan') - sd;
            
            % 3. Create Coordinates for the Shaded Area
            % Concatenate x with its flipped version, and y bounds similarly
            x_shaded = [tlin(:); flipud(tlin(:))]; % For column vectors use [x; flipud(x)] [1, 3, 6]
            y_shaded = [upper_bound(:); flipud(lower_bound(:))]; % For column vectors [upper_bound; flipud(lower_bound)]
            
            % Shade the area
            fill(x_shaded, y_shaded, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % Green, semi-transparent, no edge [3, 11, 12]
            
            
            %         plot(tlin, Fexps(:,:,jj,k), 'k:', 'linewidth', 1); hold on
            
            
            mcode = mcodes(j,:);
            
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
                modelfilename = [model,'_',cvar,'_',mvar];

                if discretized_model
                    savefolder = [modelfolder, '\Discretized\'];
                else
                    savefolder = [modelfolder, '\Approximated\'];
                end

            else
                modelfilename = [model,'_',mvar];
                savefolder = [modelfolder, '\Hill\'];
            end
                        
            for iF = iFsm
                pCa = pCai(iF,jj);
                
                if ~isnan(pCa)
                    filename = [savefolder, '\', modelfilename,'\',fibers{iF}, '\pCa=',num2str(pCa*10),'\', fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'];
                    disp(filename)
                    
                    if isfile(filename)
                        load(filename, 'tis','Cas','vis','Lis','oFi','parms', 'ts')
                    else
                        return
                    end
                    
                    oFis(iF, 1:length(oFi)) = oFi;
                end
                
                %         figure(7)
                %         plot(tis, oFi, 'color', color(j,:)); hold on
                
                
            end
            
            plot(tis - t0, mean(oFis(:,1:length(oFi)), 1, 'omitnan'), 'color', color(j,:), 'linewidth', 2); hold on
            
            if k == 1
                subplot(3,3,j+3)
                xlim([-.2 .15])
                box off
                
            else
                subplot(3,3,[7 9])
                xlim([-1.2 .15])
                box off
                
            end
            
            
        end
    end
end

%% size
figure(7)

for i = 1:3
    subplot(3,3,i)
    xlabel(gca, '     Amplitude (%L_0)      ', 'Rotation', -23, 'Position', [7 1e-5 .1])
    ylabel(gca, 'Recovery time (s)', 'Rotation', 18, 'Position', [10 1e-4 0])
end
%%

subplot(334)
text(-.3, 5.2, 'A', 'Fontsize', 16)
text(-.3, 2.2, 'B', 'Fontsize', 16)
text(-.3, -.6, 'C', 'Fontsize', 16)

subplot(332)
title('Submaximal activation: large amplitude stretch-shortening + short recovery \rightarrow reduced short-range stiffness ', 'fontsize', 8)

subplot(335)
title('Stretch-shortening + short recovery \rightarrow reduced short-range stiffness (at submaximal activation)', 'fontsize', 8)

subplot(3,3,[7 9])
title('Stretch-shortening + long recovery \rightarrow no stiffness reduction', 'fontsize', 8, 'HorizontalAlignment', 'center')
%%
figure(7)
set(gcf,'units','centimeters','position',[10 1 18 18])


end

