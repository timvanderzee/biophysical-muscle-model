% attempt to create a new time-series for figure 9

close all
iFsm = [2,3,5,6,7,8,11];
pCa = 6.1;
ISI = .001;

j = 3;

oFis = nan(max(iFsm), 2422);
mcodes = [2 2 1; 2 1 1; 1 1 1];
versions = {'parms', 'parms_v3', 'parms_v4d'};

id = [1, 2, 6];
color = acolors(id,:);

tlin = -1:.001:1;
Fexps = nan(length(tlin), max(iFsm),3);

pCas = [4.5 6.1 6.2 6.3 6.4 6.6 9];
pCai = nan(11,4);

tiso = 2;
t0 = 2*dTt + tiso;

% eth = [0 .05 .3 .7 1.2];
thmin = [.7 .3 .05 0];
thmax = [2 .7 .3 .05];

for iF = iFsm
    load([fibers{iF},'_cor_new.mat'],'data');
    
    for j = 1:length(thmin)
        
        for i = 1:length(pCas)
            pCa = pCas(i);
            [texp, Lexp, Fexp, Tsrel] = get_data(data, ISI, AMPs, pCa);
            
            Fiso = mean(Fexp(texp<-.5), 'omitnan');
            
            if Fiso < thmax(j) && Fiso > thmin(j)
                
                Fexps(:,iF,j) = interp1(texp(isfinite(Fexp)), Fexp(isfinite(Fexp)), tlin);
                pCai(iF,j) = pCa;
                break
            end
        end
    end
end

%%

pCai(:,end) = 9;
% t0 = 0;
t0 = 2*dTt + tiso;

close all
figure(1)

for j = 1:size(mcodes,1)
    subplot(1,3,j)
    for jj = 1:4
        %     subplot(1,3,jj)
        % plot(tlin, Fexps, 'k-', 'linewidth', 1); hold on
        plot(tlin, mean(Fexps(:,:,jj), 2, 'omitnan'), 'k-', 'linewidth', 2); hold on
        xlim([-.2 .15])
        box off
        
        
        %%
        
        
        mcode = mcodes(j,:);
        
        [output_mainfolder, modelfilename, ~, ~] = get_folder_and_model(mcode);
        
        for iF = iFsm
            pCa = pCai(iF,jj);
            
            if ~isnan(pCa)
                filename = [output_mainfolder{2}, '\', versions{j},'\', modelfilename,'\',fibers{iF}, '\pCa=',num2str(pCa*10),'\', fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'];
                disp(filename)
                
                load(filename, 'tis','Cas','vis','Lis','oFi','parms', 'ts')
                
                oFis(iF, :) = oFi;
            end
            
            %         figure(1)
            %         plot(tis, oFi, 'color', color(j,:)); hold on
            
            
        end
        
        plot(tis - t0, mean(oFis, 1, 'omitnan'), 'color', color(j,:), 'linewidth', 2); hold on
    end
end