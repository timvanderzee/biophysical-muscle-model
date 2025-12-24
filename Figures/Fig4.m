clear all; close all; clc
[username, githubfolder] = get_paths();
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

fig     = 6O;
savefig = 1;
iF      = 6;

acolors = get_colors;
discretized_model = 1;

% all model codes (corresponding to colors)
all_mcodes = [2 2 1; 2 1 1; 1 1 3; 1 1 2; 1 1 1; 1 2 1];

if discretized_model
    all_versions = {'', '_v3', '_v4d', '_v4d','_v4d', '_v4d'};
else
    all_versions = {'_v3', '_v3', '_v4','','_v6', '_v6'};
end

modelnames = {'Hill model (no SE)', 'Hill model (with SE)', '2-state XB model', '2-state XB coop', '3-state XB coop', '4-state XB coop'};
phases = {'Isom.','Stretch','Shortening','','Stretch','Isom.'};

%% chose figure number: specify conditions
pCas = [4.5 6.2 9; 4.5 6.2 9; 4.5 6.2 9];

% chosen ISIs, AMPs and pCas
if fig == 4 || fig == 5
    ISIs = [.1 .1 .1; .001 .001 .001; .1 .1 .1;]; 
    AMPs = [.0383 .0383 .0383; 0 0 0; .0383 .0383 .0383];
    tmaxs =  [0.4; .14; 0.4];
    zylim = [30 60; 15 45];
    jid = 1;
    
elseif fig == 6
    ISIs = [ .001 .001 .001; .316 .316 .316; .001 .001 .001];
    AMPs = [.0383 .0383 .0383; .0383 .0383 .0383; .0383 .0383 .0383]; 
    tmaxs = [0.32 0.6345 0.32];
    zylim = [30 60; 5 35];
    jid = 1;
end

dTt = .0383/.4545; % test stretch (= constant)
dTc = .0383 / .4545; % conditioning stretch
tiso = dTt*3+dTc*2+ISIs(1) + 2;
Ts = [tiso dTc dTc ISIs(1)];
dt = tiso - sum(Ts);

zxlim = [-.01 .051; -.01 .051] + [0; -dt];
ls = {'-',':','-'};

%% choose fiber: load data and parameters
% chosen models
if fig == 4
    mid = [1 2 3];
    
elseif fig == 5
    mid = [3 4 5];

elseif fig == 6
    mid = [4 5 6];
end

colors      = acolors(mid,:);
versions    = all_versions(mid);
mcodes      = all_mcodes(mid,:);

%% load model

for j = 1:size(ISIs,1)
    for i = 1:size(ISIs,2)

        ISI = ISIs(j,i);
        AMP = AMPs(j,i);

        for kk = 1:size(mcodes,1)
            mcode = mcodes(kk,:);
            [output_mainfolder, modelfilename, ~, ~] = get_folder_and_model(mcodes(kk,:));
            
            filename = [output_mainfolder{2}, '\parms', versions{kk},'\', modelfilename,'\',fibers{iF}, '\pCa=',num2str(pCas(j,i)*10),'\', fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'];
            disp(filename)

            if ~exist(filename, 'file')
                disp('using approximated')
                filename = [output_mainfolder{2}, '\parms_v4\', modelfilename,'\',fibers{iF}, '\pCa=',num2str(pCas(j,i)*10),'\', fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'];
            end

            Mdata(j,i,kk) = load(filename, 'tis','Cas','vis','Lis','oFi','parms', 'ts');
        end
    end
end

%% make axes
figure(1)
set(gcf,'units','centimeters','position',[10 2 20 15])

tmax = max(tmaxs);

for i = 1:8
    ax(i) = axes('Xlim', [-.05, tmax], 'Xtick', (-.05:.05:tmax), 'box', 'off', 'Fontsize', 6);
    hold(ax(i), 'on');
    
    if i > 1
        set(ax(i), 'Ytick', 0:10:200)
    end
end

ylabel(ax(1), 'Length (%L_0)', 'Fontsize', 8)
set(ax(1), 'Ylim', [-.5 5])

xlabel(ax(2), 'Time (s)', 'Fontsize', 8)
xlabel(ax(5), 'Time (s)', 'Fontsize', 8)
xlabel(ax(8), 'Time (s)', 'Fontsize', 8)
ylabel(ax(2), 'Force (%F_0)', 'Fontsize', 8)
set(ax(2), 'Ylim', [0 200])

for kk = 1:2
    for k = 1:3
        set(ax(k+2+(kk-1)*3), 'xlim', zxlim(kk,:), 'ylim', zylim(kk,:))
        
        if (k+2+(kk-1)*3) < 6
            ylabel(ax(k+2+(kk-1)*3), 'Force (%F_0)', 'fontsize', 8)
        end
    end
end

stretches = {'Cond.', 'Test'};
title(ax(3), stretches{1}, 'fontsize', 8)
title(ax(6), stretches{2}, 'fontsize', 8)

% width should depend on time range
fac = 14;

set(ax(1), 'units','centimeters', 'Position', [2 12.5 tmax*fac 1.5], 'xticklabels', {})
set(ax(2), 'units','centimeters', 'Position', [2 2 tmax*fac 10])

% left-most insets
set(ax(3), 'units','centimeters', 'Position', [tmax*fac + 3 9 1.5 3], 'xticklabels', {});
set(ax(4), 'units','centimeters', 'Position', [tmax*fac + 3 5.5 1.5 3], 'xticklabels', {});
set(ax(5), 'units','centimeters', 'Position', [tmax*fac + 3 2 1.5 3]);

% right-most insets
set(ax(6), 'units','centimeters', 'Position', [tmax*fac + 5 9 1.5 3], 'xticklabels', {});
set(ax(7), 'units','centimeters', 'Position', [tmax*fac + 5 5.5 1.5 3], 'xticklabels', {});
set(ax(8), 'units','centimeters', 'Position', [tmax*fac + 5 2 1.5 3]);

for kk = 1:2
    for k = 1:3
        rectangle(ax(k+2+(kk-1)*3), 'Position', [zxlim(kk,1) zylim(kk,1) diff(zxlim(kk,:)) diff(zylim(kk,:))])
    end

    text(ax(2), mean(zxlim(kk,:)), zylim(kk,2)+5, stretches{kk}, 'fontsize', 6, 'HorizontalAlignment', 'center', 'fontweight', 'bold');
    rectangle(ax(2), 'Position', [zxlim(kk,1) zylim(kk,1) diff(zxlim(kk,:)) diff(zylim(kk,:))])
end

%% plot model
% brigthen factor
bf = [0 .5 0];

for j = 1:size(ISIs,1)
    for i = 1:size(ISIs,2) 
        
        dTt = .0383/.4545; % test stretch (= constant)
        dTc = AMPs(j,i) / .4545; % conditioning stretch
        tiso = dTt*3+dTc*2+ISIs(j,i) + 2;
        Ts = [tiso dTc dTc ISIs(j,i)];
        dt = tiso - sum(Ts);
        
        for kk = 1:size(mcodes,1)
        
            t = Mdata(j,i,kk).tis + 3*dTt - tiso;
           
            t = t - dt;
            tmax2 = .15 - dt;

            plot(ax(2), t(t<tmax2), Mdata(j,i,kk).oFi(t<tmax2)*100, 'linestyle', ls{j}, 'linewidth',2, 'color', brighten(colors(kk,:), bf(j))); hold on
            
            if j == 1
                for k = 1:2
                    plot(ax(kk+2 + (k-1)*3), t(t<tmax2), Mdata(j,i,kk).oFi(t<tmax2)*100, 'linestyle', ls{j}, 'linewidth',2, 'color', brighten(colors(kk,:), bf(j))); hold on
                end
            end
        end       
    end
end


%% plot data
% load data
[output_mainfolder, filename, ~, ~] = get_folder_and_model(mcodes(1,:));
cd([output_mainfolder{2},'\data'])
load([fibers{iF},'_cor_new.mat'],'data');

for j = 1:size(ISIs,1)

    % extra data
    [texp, Lexp, Fexp, Tsrel] = get_data(data, ISIs(j,:), AMPs(j,:), pCas(j,:));
    
    dt(j) = Tsrel(1,2);
    
    texp = texp - dt(j) - .003;
        
    % add vertical lines
    if j == jid
        tids = sort([Tsrel(1,:) .15]) - Tsrel(1,2);

        for ii = 1:(length(tids)-2)
            for k = 1:8
                xline(ax(k), tids(ii), ':', 'color', [.5 .5 .5]); hold on
            end

            text(ax(1), mean([tids(ii) tids(ii+1)]), 4.5, phases{ii}, 'Fontsize',7,'HorizontalAlignment','center')
        end
    end

    % add data
    for i = 1:size(ISIs,2)
    
        % only certain samples
        id = texp(:,i) < tmaxs(j);
        
        tmax = max(tmaxs);
        
        % length
        plot(ax(1), texp(id,i), Lexp(id,i)*100, 'color', brighten([.5 .5 .5], bf(j)), 'linewidth', 2, 'linestyle',ls{j}); 
        plot(ax(2), texp(id,i), Fexp(id,i)*100, 'color', brighten([.5 .5 .5], bf(j)), 'linewidth', 2, 'linestyle',ls{j}); 
        
        % extra panels
        if j == 1
            for k = 1:2
                for kk = 1:3
                    plot(ax(kk+2 + (k-1)*3), texp(id,i), Fexp(id,i)*100, 'color', brighten([.5 .5 .5], bf(j)), 'linewidth', 2, 'linestyle',ls{j}); 
                end
            end
        end
    end
end

% legend
h = legend(ax(2), 'labels', modelnames(mid),'location','bestoutside', 'Fontsize', 8);
set(h, 'units','centimeters', 'Position', [tmax*fac + 3 13 3.5 1])

%% optionally export to PNG
cd(['C:\Users\',username,'\OneDrive\9. Short-range stiffness\figures\MAT'])

if discretized_model
    figname = ['Fig',num2str(fig),'.png'];
else
    figname = ['FigS',num2str(fig-3),'.png'];
end


if savefig
    figure(1)
    exportgraphics(gcf,figname)
end







