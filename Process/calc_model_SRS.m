clear all; close all; clc
[username, githubfolder] = get_paths();

mcodes = [2 1 1; 1 1 1; 1 1 3; 1 2 1];
mcodes = [1 1 1];
iFs = [1 2 3, 5,6, 7, 8, 10, 11];
AMPs = [0 12 38 121 216 288 383 532 682]/10000;
ISIs = [1 10 50 100 200 316 500 1000 3160 10000]/1000;
pCas = [4.5 6.1 6.2 6.3 6.4 6.6 9];
Ca = 10.^(-pCas+6);
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

% visualize = 0;

% AMPs = 532 / 10000;

for iii = 1:size(mcodes,1)

% load parameters
mcode = mcodes(iii,:);

% load parameters
% mcode = [1 1 3];
[output_mainfolder, filename, ~, ~] = get_folder_and_model(mcode);

cd([githubfolder, '\biophysical-muscle-model\Parameters'])
load(['parms_',filename,'.mat'], 'pparms')


%% step 2: stretch-shortening
% conditions
% AMPs = [0 12 38 121 216 288 383 532 682]/10000;
% % ISIs = [1 10 100 316 1000]/1000;
% ISIs = [1 10 100 316 1000 3160 10000]/1000;
% pCas = [4.5 6.1 6.2 6.3 6.4 6.6 9];
% Ca = 10.^(-pCas+6);

% outputs
F0      = nan(length(pCas), length(ISIs),length(AMPs), 11);
Scond   = nan(length(pCas), length(ISIs),length(AMPs), 11);
Stest   = nan(length(pCas), length(ISIs),length(AMPs), 11);

visualize = 0;

for iF = iFs
    for i = 1:length(Ca)
        for ii = 1:length(AMPs)
            
            AMP = AMPs(ii);
            dTt = .0383/.4545; % test stretch (= constant)
            dTc = AMP / .4545; % conditioning stretch
            
            for jj = 1:length(ISIs)
                
                ISI = ISIs(jj);
                
                tiso = dTt*3+dTc*2+ISI + 2;
                
                cd([output_mainfolder{2}])
                
                cd([filename,'\',fibers{iF}, '\pCa=',num2str(pCas(i)*10)])
                
                disp([fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'])
                load([fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'], ...
                    'tis','Cas','vis','Lis','oFi','parms', 'ts')
                
                [id0,id1,id2] = get_indices(tis, tiso, ts, dTt, dTc, ISI, Ca(i));
                
                if visualize
                    subplot(411)
                    plot(tis, Cas, 'b', 'linewidth',1);
                    
                    subplot(412)
                    plot(tis, vis,'b',  'linewidth',1);
                    
                    subplot(413)
                    plot(tis, Lis,'b',  'linewidth',1);
                    
                    subplot(414);
                    plot(tis, oFi,'b',tis(id1), oFi(id1),'k.')
                    
                    if AMP == .0383
                        hold on
                        plot(tis(id2), oFi(id2),'g.', tis(id0), oFi(id0),'r.')
                        hold off
                    end
                    
                    for j = 1:4
                        subplot(4,1,j)
                        box off
                    end
                end
                
                np1 = polyfit(Lis(id1), oFi(id1), 1);
                Stest(i,jj,ii,iF) = np1(1);
                
                if AMP == .0383
                    np2 = polyfit(Lis(id2), oFi(id2), 1);
                    Scond(i,jj,ii,iF) = np2(1);
                    F0(i,jj,ii,iF) = mean(oFi(id0));
                end
                
            end
        end
    end
end




%% save
cd(githubfolder)
cd('biophysical-muscle-model')
cd('Model output\SRS')
save([filename, '_SRS.mat'])


end

return
%% plot for all fibers
close all
SRSrel = Stest./Scond(:,:,AMPs == .0383,:);

aAMPs = repmat(AMPs, 7, 1);
aISIs = repmat(ISIs(:), 1, 8);

iF = 6;
for i = 1:size(SRSrel,1)
    SRS = squeeze(SRSrel(i,:,:,iF));
    
    figure(1)
    nexttile
    surf(aAMPs, aISIs, SRS)
    set(gca, 'Yscale', 'log')
    
end

