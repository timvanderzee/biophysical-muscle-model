clear all; close all; clc
[username, githubfolder] = get_paths();

Kss = 1:7;
tiso = 3;

% iFs = [1 2 3, 5, 6, 7, 8, 10, 11];
iFs = [2,3,5,6,7,8,11];
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
mcodes = [2 1 1; 1 1 1; 1 1 3; 1 2 1];
mcodes = [1 1 1];

visualize = 0;

%% calc RMSD
AMPs = [0 12 38 121 216 288 383 682]/10000;
ISIs = [1 10 100 316 1000 3160 10000]/1000;
pCas = [4.5 6.1 6.2 6.3 6.4 6.6 9];
Ca = 10.^(-pCas+6);

for iii = 1:size(mcodes,1)
    RMSD = nan(length(pCas), length(ISIs), length(AMPs), iFs(end), 7);
    
    mcode = mcodes(iii,:);
    [output_mainfolder, modelname, ~, ~] = get_folder_and_model(mcode);
    
    for iF = iFs
        cd(['C:\Users\',username,'\OneDrive - KU Leuven\9. Short-range stiffness\matlab\data'])
        load([fibers{iF},'_cor_new.mat'],'data')
        
        for i = 1:length(Ca)
            %         cd([output_mainfolder{2}])
            cd([output_mainfolder{2}, '\parms_v4'])
            
            cd([modelname,'\',fibers{iF}, '\pCa=',num2str(pCas(i)*10)])
            
            for m = 1:length(AMPs)
                
                AMP = AMPs(m);
                dTt = .0383/.4545; % test stretch (= constant)
                dTc = AMP / .4545; % conditioning stretch
                
                for n = 1:length(ISIs)
                    ISI = ISIs(n);
                    
                    filename = [fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'];
                    disp(filename)
                    
                    if exist(filename, 'file')
                        
                        try
                            load(filename, 'tis','Cas','vis','Lis','oFi','parms', 'ts')
                            
                            tiso = dTt*3+dTc*2+ISI;
                            
                            texp = data.texp(:,i,n,m) - .005;
                            Fexp = data.Fexp(:,i,n,m);
                            
                            % center around second stretch
                            tm = tis - 2 -ISI - 2 * dTc;
                            oFii = interp1(tm, oFi, texp);
                            
                            tids = [-ISI - 2*dTc - .1; -ISI - 2*dTc; -ISI - dTc; -ISI; 0; dTt; .16];
                            
                            for ii = 1:length(tids)
                                if ii < length(tids)
                                    id = texp < tids(ii+1) & texp >= tids(ii);
                                else % overall
                                    id = texp < tids(end) & texp >= tids(1);
                                end
                                
                                if sum(id) > 0
                                    % compute RMSD
                                    %                 id = tm < .16 & tm > (-ISI - 2 * dTc - .1);
                                    if visualize
                                        figure(1)
                                        plot(tm, oFi, '-', texp, Fexp, '-', texp(id), oFii(id), '.')
                                        pause
                                    end
                                    
                                    RMSDs = sqrt((oFii(id) - Fexp(id)).^2) * 100;
                                    
                                    RMSD(i,n,m,iF,ii) = sqrt(mean((oFii(id) - Fexp(id)).^2, 'omitnan'));
                                else
                                    RMSD(i,n,m,iF,ii) = nan;
                                end
                            end
                            
                        catch
                            disp(['Unable to read: ', filename])
                        end
                        
                    end
                    
                end
            end
        end
    end
    
    
    %% save
    cd(githubfolder)
    cd('biophysical-muscle-model/Model output/RMSD')
    
    save([modelname, '_RMSD.mat'])
end

return

%%
close all
aAMPs = repmat(AMPs, 7, 1);
aISIs = repmat(ISIs(:), 1, 8);

iF = 2;
for i = 1:size(RMSD,1)
    eps = squeeze(RMSD(i,:,:,iF));
    
    figure(1)
    nexttile
    surf(aAMPs, aISIs, eps)
    set(gca, 'Yscale', 'log')
    
end



return
%%

RMSDm = mean(RMSD, 5, 'omitnan');
RMSDi = nan(size(RMSDm));
RMSDc = nan(size(RMSDm));



% correct for individual offset
for iF = 1:9
    RMSDc(:,:,iF,:) = RMSDm(:,:,iF,:) - mean(RMSDm(:,:,iF,:),'all','omitnan') + mean(RMSDm(:),'all','omitnan');
end

return

%%

close all
figure(1)

% for i = 1:size(AMPs,3)
%     nexttile
i = 1;

surf(AMPs(:,:,i), ISIs(:,:,i), mean(RMSDc,3,'omitnan'));

set(gca,'YScale','log')
xlabel('Amplitude')
ylabel('Recovery time')

% title(['pCa = ', num2str(mean(ACTs(:,:,i),'all'))])

hold on

plot3(AMPs(1,1,i), ISIs(1,1,i), mean(RMSDc(1,1,:),3,'omitnan'),'ro')
plot3(AMPs(4,2,i), ISIs(4,2,i), mean(RMSDc(4,2,:),3,'omitnan'),'ro')

view(45, 45)
zlim([0 .05])
% end

%%

close all

ms = [4 4 2 2];
ns = [1 4 4 1];

titles = {'Overall','Short recovery','Long recovery','Long recovery','Short recovery'};
subtitles = {'','Large amplitude','Large amplitude','Small amplitude','Small amplitude'};

mnames = {'Hill','No coop.','Coop.', 'Coop. + FD'};

pcolors = flip(parula(7));

figure(1)
color = get(gca,'colororder');
colors = [color(2,:); pcolors(4:end-1,:)];

jid = 1;
jr = 1;

X = nan(length(jid), size(RMSDc,3));
Y = nan(length(jid), size(RMSDc,3));

mm = 1;

for i = 1:4
    subplot(1,5,i+1)
    
    %     v = violinplot(squeeze(RMSDi(ms(i),ns(i),:,jid) - RMSDi(ms(i),ns(i),:,jr))*100, mnames, 'ViolinColor', colors,'Edgecolor', [1 1 1]);
    v = violinplot(squeeze(RMSDc(ms(i),ns(i),:,jid))*100, mnames, 'ViolinColor', colors,'Edgecolor', [1 1 1],'ShowMean',true, 'ShowMedian',false);
    
    for j = 1:length(jid)
        X(j,:) = v(j).ScatterPlot.XData;
        Y(j,:) = v(j).ScatterPlot.YData;
    end
    
    plot(X, Y, ':','color',[.5 .5 .5])
    
end

%%
subplot(151)
% v = violinplot((mmRMSDi(:,jid)-mmRMSDi(:,jr))*100, mnames, 'ViolinColor', colors,'Edgecolor', [1 1 1]);
v = violinplot((mmRMSDc(:,jid))*100, mnames, 'ViolinColor', colors,'Edgecolor', [1 1 1],'ShowMean',true, 'ShowMedian',false);

for j = 1:length(jid)
    X(j,:) = v(j).ScatterPlot.XData;
    Y(j,:) = v(j).ScatterPlot.YData;
end

plot(X, Y, ':','color',[.5 .5 .5])

drawnow

%
figure(mm)

for i = 1:5
    subplot(1,5,i)
    %     axis([.5 4.5 -100 100])
    axis([.5 4.5 -90 60])
    %         axis([.5 3.5 -4 4])
    box off
    %     axis tight
    title(titles{i})
    subtitle(subtitles{i},'Fontweight','bold')
    ylabel('\DeltaRMSD (%)')
    
    yline(0,'k-','linewidth',1.5)
end

subplot(151)
% a = annotation('textarrow',[.21 .21],[.7 .52],'String','Hill-type model','Fontsize',8,'HeadWidth',5,'HeadLength',5);

% make nice
for i = 2:5
    subplot(1,5,i)
    ylabel('')
    set(gca,'YTickLabel', '')
    
end
%
figure(mm)
set(gcf,'units','normalized','position',[.2 .2 .6 .3])

ABC = {'A','B','C','D','E'};
for i = 1:5
    subplot(1,5,i)
    ym = max(get(gca,'ylim'));
    text(0,ym*1.2, ABC{i},'fontsize',16)
end


%%
close all
figure(1)

for kk = 1:size(mcodes,1)
    nexttile
    surf(ISIs, AMPs, RMSD(:,:,kk))
    xlabel('ISI')
    ylabel('AMP')
    
    set(gca, 'XScale', 'log')
    % zlim([0 2])
end