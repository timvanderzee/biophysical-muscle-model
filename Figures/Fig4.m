clear all; close all; clc
[username, githubfolder] = get_paths();
savefig = 0;

acolors = lines(6);
acolors(1,:) = brighten(acolors(1,:), .2);
acolors(2,:) = brighten(acolors(2,:), -.7);
acolors(3,:) = brighten(acolors(3,:), .5);
acolors(4,:) = brighten(acolors(4,:), .5);
acolors(5,:) = brighten(acolors(5,:), -.5);
acolors(6,:) = brighten(acolors(6,:), .7);

discretized_model = 1;

%% chose figure number: specify conditions
fig = 5;
iF = 6;

pCas = [4.5 6.2 9;
        4.5 6.2 9;
        4.5 6.2 9];

% chosen ISIs, AMPs and pCas
if fig == 4 || fig == 5
    ISIs = [.1 .1 .1;
            .001 .001 .001;
            .1 .1 .1;]; 
    
    AMPs = [.0383 .0383 .0383;
            0      0 0;
            .0383 .0383 .0383];% solid
    
    
    titles = {'Maximal activation', 'Submaximal activation'};
    
    tmaxs =  [0.4; .14; 0.4];
    
elseif fig == 6
    
    ISIs = [ .316 .316 .316;
            .001 .001 .001;
             .316 .316 .316];
    
    AMPs = [.0383 .0383 .0383;
            .0383 .0383 .0383;
            .0383 .0383 .0383];
        
    tmaxs = [0.6345 0.32 0.6345];
    
    titles = {'Maximal activation', 'Submaximal activation'};
    
end


%% choose fiber: load data and parameters
if fig == 4
    mcodes = [2 2 1; 2 1 1; 1 1 3];
    colors = [acolors(end,:); acolors];
    
    %     versions = {'_v3', '_v3', '_v4'};
    
    if discretized_model
        versions = {'', '_v3', '_v4d', '_v4d'};
    else
        versions = {'_v3', '_v3', '_v4'};
    end
    
elseif fig == 5
    mcodes = [1 1 3; 1 1 2; 1 1 1;  1 2 1];
    colors = acolors(2:end,:);
   
    
    if discretized_model
        versions = {'_v4d', '_v4d', '_v4d', '_v4d'};
    else
        versions = {'_v3', '_v3', '_v4'};
    end
    
else
    mcodes = [1 1 3; 1 1 2; 1 1 1; 1 2 1];
    colors = acolors(2:end,:);
    
    
    if discretized_model
        versions = {'_v3', '_v4d', '_v4d','_v4d'};
    else
%         versions = {'_v4', '_v4'};
        versions = {'_v4', '_v6'};
    end
end

% load data
[output_mainfolder, filename, ~, ~] = get_folder_and_model(mcodes(1,:));
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
cd([output_mainfolder{2},'\data'])
load([fibers{iF},'_cor_new.mat'],'data');

%% make axes
figure(1)

for i = 1:5
    ax(i) = axes();
end


%% evaluate
odeopt = odeset('maxstep', 1e-2);
ls = {'-',':','-'; '-',':','-'; '-',':','-'; '--',':','--'};

bf = [0 .5 0];

for j = 1:size(ISIs,1)

    [texp, Lexp, Fexp, Tsrel] = get_data(data, ISIs(j,:), AMPs(j,:), pCas(j,:));
    texp = texp - Tsrel(1,2);

    for i = 1:size(ISIs,2)
    
        id = texp(:,i) < tmaxs(j);
        
        tmax = max(tmaxs);
        

        plot(ax(1), texp(id,i), Lexp(id,i)*100, 'color', brighten([.5 .5 .5], bf(j)), 'linewidth', 2, 'linestyle',ls{1,j}); 
        hold(ax(1), 'on')
        axis(ax(1), [-.05 tmax -.5 5])
        set(ax(1), 'box', 'off', 'Fontsize', 6)
        xticks(ax(1), -.05:.05:tmax);
        ylabel(ax(1), 'Length (%L_0)', 'Fontsize', 8)
        
        for k = 1:4
            plot(ax(k+1), texp(id,i), Fexp(id,i)*100, 'color', brighten([.5 .5 .5], bf(j)), 'linewidth', 2, 'linestyle',ls{1,j}); 
            hold(ax(k+1), 'on')
            axis(ax(k+1), [-.05 tmax 0 200])
            xticks(ax(k+1), -.05:.05:tmax);

            set(ax(k+1), 'box', 'off', 'Fontsize', 6)
        end
        
        ylabel(ax(2), 'Force (%F_0)', 'Fontsize', 8)
        
        
        % evaluate fit
        Ca = 10.^(-pCas(j,i)+6);
        
        dTt = .0383/.4545; % test stretch (= constant)
        dTc = AMPs(j,i) / .4545; % conditioning stretch
        ISI = ISIs(j,i);
        AMP = AMPs(j,i);
        
        tids = sort([Tsrel(i,:) .15]) - Tsrel(1,2);
        
        for kk = 1:size(mcodes,1)
            mcode = mcodes(kk,:);
            [output_mainfolder, modelfilename, ~, ~] = get_folder_and_model(mcodes(kk,:));
            tiso = dTt*3+dTc*2+ISI + 2;
            
            filename = [output_mainfolder{2}, '\parms', versions{kk},'\', modelfilename,'\',fibers{iF}, '\pCa=',num2str(pCas(j,i)*10),'\', fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'];
            disp(filename)
            
            if ~exist(filename, 'file')
                disp('using approximated')
                filename = [output_mainfolder{2}, '\parms_v4\', modelfilename,'\',fibers{iF}, '\pCa=',num2str(pCas(j,i)*10),'\', fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'];
                
            end
            
            load(filename, 'tis','Cas','vis','Lis','oFi','parms', 'ts')
            
            t = tis + 3*dTt - tiso;
            
            t = t - Tsrel(1,2);
            tmax2 = .15 - Tsrel(1,2);
            
            for k = 1:4
                plot(ax(k+1), t(t<tmax2), oFi(t<tmax2)*100, 'linestyle', ls{kk,j}, 'linewidth',2, 'color', brighten(colors(kk,:), bf(j))); hold on
            end
            
            % compute RMSD
            id = t < .15 & t > (-ISI - 2 * dTc - .1);
            oFii = interp1(t(id), oFi(id), texp(:,i));
            
            xlabel('Time (s)', 'Fontsize', 8)

        end
        
        % add vertical lines
        if ISI > .01
            phases = {'Isometric','Pre-stretch','Shortening','Recovery','Test stretch','Isom.'};
        else
            phases = {'Isometric','Pre-stretch','Shortening','','Test stretch','Isom.'};
        end
        

        for ii = 1:(length(tids)-2)
            xline(ax(1), tids(ii), ':', 'color', [.5 .5 .5]); hold on
            
            if j == 3 && i == 1
                text(ax(1), mean([tids(ii) tids(ii+1)]), 4.5, phases{ii}, 'Fontsize',7,'HorizontalAlignment','center')
            end
        end
        
%         ax2 = subplot(5,1,[2 5])
        for ii = 1:(length(tids)-2)
            xline(ax(2), [tids(ii) tids(ii)], ':', 'color', [.5 .5 .5]); hold on
        end
        
        
    end
end

%%
figure(1)

if fig == 4
    h = legend(ax(2), 'Data','Hill (no SE)','Hill (with SE)', '2-state','location','bestoutside', 'Fontsize', 8);
elseif fig == 5
    h = legend(ax(2), 'Data','2-state','2-state coop','3-state coop','4-state coop','location','bestoutside', 'Fontsize', 8);
else
    h = legend(ax(2), 'Data','2-state','2-state coop','3-state coop','4-state coop','location','bestoutside', 'Fontsize', 8);
end

trange = [-.05 tmax];
set(gcf,'units','centimeters','position',[10 2 15 15])

%%
for k = 1:3
    set(ax(k+2), 'xlim', [.25 .35], 'ylim', [20 60])
end

%% size
figure(1)

set(ax(1), 'units','centimeters', 'Position', [2 12 8 2])
set(ax(2), 'units','centimeters', 'Position', [2 2 8 9])
set(h, 'units','centimeters', 'Position', [11 12 3 2])

% ax3 = axes();
% ax4 = axes();
% ax5 = axes();

%% 
set(ax(3), 'units','centimeters', 'Position', [11 2 3 2.5], 'xticklabels', {}, 'yticklabels', {});
set(ax(4), 'units','centimeters', 'Position', [11 5 3 2.5], 'xticklabels', {}, 'yticklabels', {});
set(ax(5), 'units','centimeters', 'Position', [11 8 3 2.5], 'xticklabels', {}, 'yticklabels', {});

%%
plot(ax(3), texp(:,i), Fexp(:,i)*100, 'color', brighten([.5 .5 .5], bf(j)), 'linewidth', 2, 'linestyle',ls{1,j}); 

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







