function[] = normalize_data(force_pCa_folder, reorganized_datafolder, outputfolder)

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

cd(force_pCa_folder)
F_pCa = load('force_pCa.mat','F0','Fmax','pCas');
pCas = F_pCa.pCas;
extra_plots = 0;

if extra_plots
    %% one plot per pCa
    if ishandle(1), close(1); end
    
    cd(reorganized_datafolder)
    
    colors = parula(length(fibers));
    for k = 1:length(fibers)
        load([fibers{k},'_reorganized.mat'],'Time','Force','pCas', 'Length')
        
        disp(fibers{k})
        figure(1)
        
        for i = 1:length(pCas)
            subplot(4,2,i)
            plot(Time(:,i,3,7),Force(:,i,3,7),'color',colors(k,:),'linewidth',1.5); hold on
            title(pCas(i))
            if isfinite(F_pCa.F0(k, F_pCa.pCas == round(pCas(i),1)))
                yline(F_pCa.F0(k, F_pCa.pCas == round(pCas(i),1)),'--','color',colors(k,:))
            end
        end
        
    end
    
    %% one plot per pCa, normalized
    
    if ishandle(2), close(2); end
    figure(2)
    
    colors = parula(length(fibers));
    for k = 1:length(fibers)
        load([fibers{k},'_reorganized.mat'],'Time','Force','pCas')
        
        disp(fibers{k})
        figure(2)
        
        for i = 1:length(pCas)
            subplot(4,2,i)
            plot(Time(:,i,3,7),Force(:,i,3,7)/F_pCa.Fmax(k),'color',colors(k,:),'linewidth',1.5); hold on
            title(pCas(i))
            ylim([0 2])
            
            if isfinite(F_pCa.F0(k, F_pCa.pCas == round(pCas(i),1)))
                yline(F_pCa.F0(k, F_pCa.pCas == round(pCas(i),1))/F_pCa.Fmax(k),'--','color',colors(k,:))
            end
        end
    end
    
    %% one plot per fiber
    if ishandle(3), close(3); end
    figure(3)
    
    colors = parula(length(pCas));
    
    for k = 1:length(fibers)
        
        load([fibers{k},'_reorganized.mat'],'Time','Force','pCas')
        nexttile
        
        for i = 1:length(pCas)
            
            plot(Time(:,i,3,7),Force(:,i,3,7),'color',colors(i,:),'linewidth',1.5); hold on
            ylim([0 500])
            box off
            xlabel('Time (s)')
            ylabel('Force (kPa)')
            title(fibers{k})
            xlim([-.5 .5])
            
            plot(-.5, F_pCa.F0(k, F_pCa.pCas == round(pCas(i),1)),'o','color',colors(i,:),'markerfacecolor',colors(i,:))
            
        end
    end
    
    %% one plot per fiber, normalized
    colors = parula(length(pCas));
    
    if ishandle(4), close(4); end
    figure(4)
    
    for k = 1:length(fibers)
        
        load([fibers{k},'_reorganized.mat'],'Time','Force','pCas')
        nexttile
        
        for i = 1:length(pCas)
            
            plot(Time(:,i,3,7),Force(:,i,3,7)/F_pCa.Fmax(k),'color',colors(i,:),'linewidth',1.5); hold on
            ylim([0 2])
            box off
            xlabel('Time (s)')
            ylabel('Force (-)')
            title(fibers{k})
            xlim([-.5 .5])
            
            plot(-.5, F_pCa.F0(k, F_pCa.pCas == round(pCas(i),1))/F_pCa.Fmax(k),'o','color',colors(i,:),'markerfacecolor',colors(i,:))
        end
    end
end


%% one plot per fiber, corrected
colors = parula(length(pCas));
cd(reorganized_datafolder)
    
if ishandle(4), close(4); end
figure(4)

load([fibers{1},'_reorganized.mat'],'Force')
Force_cor   = nan(size(Force,1), size(Force,2), length(fibers), size(Force,3), size(Force,4));
Time_cor    = nan(size(Force,1), size(Force,2), length(fibers), size(Force,3), size(Force,4));
Length_cor    = nan(size(Force,1), size(Force,2), length(fibers), size(Force,3), size(Force,4));

fac = nan(size(Force,2), length(fibers), size(Force,3), size(Force,4));

AMPs = [0 12 38 121 216 288 383 682];
ISIs = [1 10 100 316 1000 3162 10000];

show = 0;
F_pCa.F02 =  nan(size(Force,2), length(fibers), size(Force,3), size(Force,4));
 
for k = 1:length(fibers)
    disp([fibers{k},'_reorganized.mat'])
    load([fibers{k},'_reorganized.mat'],'Time','Force','pCas','Length')
    
    if show
        nexttile
    end
    
    for n = 1:length(ISIs) % 3
        for m = 1:length(AMPs) % 7
            
            % remove NaNs (replicable the pCa = 9)
            Time_cor(:,:,k,n,m) = repmat(Time(:,end,n,m),1,length(pCas));
            
            for i = 1:length(pCas)
                
                % before the first stretch
                id = Time(:,i,n,m) < -((ISIs(n) + 2*AMPs(m)/4.545)/1000);
                
                % normalize the length
                Length_cor(:,i,k,n,m) = (Length(:,i,n,m)-mean(Length(id,i,n,m),'omitnan'))/mean(Length(id,i,n,m),'omitnan');
                
                % correct the force
                Force_cor(:,i,k,n,m) = Force(:,i,n,m)/mean(Force(id,i,n,m),'omitnan') * F_pCa.F0(k,i)/F_pCa.Fmax(k);
                
                F_pCa.F02(i,k,n,m) = mean(Force(id,i,n,m),'omitnan');
                fac(i,k,n,m) = mean(Force(id,i,n,m),'omitnan') / F_pCa.F0(k,i);
                
                if show
                    % plot the corrected force
                    plot(Time_cor(:,i,k,n,m), Force_cor(:,i,k,n,m),'color',colors(i,:),'linewidth',1.5); hold on
                    
                    ylim([0 2])
                    box off
                    xlabel('Time (s)')
                    ylabel('Force (-)')
                    title(fibers{k})
                    %                 xlim([-.5 .5])
                    
                    plot([Time(1,i,n,m) -.5], ones(1,2) * F_pCa.F0(k, F_pCa.pCas == round(pCas(i),1))/F_pCa.Fmax(k),'o','color',colors(i,:),'markerfacecolor',colors(i,:))
                end
            end
        end
    end
end

%% plot all fibers in one figure
colors = parula(length(fibers));


n = 3;
m = 7;


if ishandle(5), close(5); end
figure(5)

for i = 1:length(pCas)
    nexttile
    
    for k = 1:length(fibers)
        plot(Time_cor(:,i,k,n,m), Force_cor(:,i,k,n,m), 'color', colors(k,:), 'linewidth', 1.5); hold on
        title(num2str(pCas(i)));
        box off
    end
end

%% plot force pCa
% figure(6)
% plot(F_pCa.pCas, F_pCa.F0./F_pCa.Fmax,'o')


%% save fit data: 0.1 s, 0.038 L0
% n = 3;
% m = 7;
%
% ns = [3 3];
% ms = [1 7];
%
% close all
%
% cd('C:\Users\timvd\Documents\muscle-thixotropy\data')
%
% for i = 1:2
%     n = ns(i);
%     m = ms(i);
%
%     data(i).texp = Time_cor(:,:,:,n,m);
%     data(i).Fexp = Force_cor(:,:,:,n,m);
%     data(i).pCas = pCas;
%     data(i).Act = F_pCa.F0./F_pCa.Fmax;
%
%     data(i).AMP = [AMPs(m)/10000 .0383];
%     data(i).ISI = ISIs(n) / 1000;
%
% end
% save('All_fibers_fit_data.mat','data')

%% save all data (too big for GitHub)
cd(outputfolder)

for k = 1:length(fibers)
    
    data.texp = squeeze(Time_cor(:,:,k,:,:));
    data.Fexp = squeeze(Force_cor(:,:,k,:,:));
    data.pCas = pCas;
    data.AMPs = AMPs;
    data.ISIs = ISIs;
    data.Act = F_pCa.F0(k,:)./F_pCa.Fmax(k);
    
    data.Lexp = squeeze(Length_cor(:,:,k,:,:));
    
    save([fibers{k},'.mat'],'data')
    
end


