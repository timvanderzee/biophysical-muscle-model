clear all; close all; clc
[username, githubfolder] = get_paths();

Kss = 1:7;
tiso = 3;

% iFs = [1 2 3, 5, 6, 7, 8, 10, 11];
iFs = 5 %[2,3,5,6,7,8,11];
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
mcodes = [2 1 1; 1 1 1; 1 1 3; 1 2 1];
mcodes = [1 1 1];

visualize = 0;

AMPs = [0 12 38 121 216 288 383 682]/10000;
ISIs = [1 10 100 316 1000 3160 10000]/1000;
pCas = [4.5 6.1 6.2 6.3 6.4 6.6 9];
% 
% pCas = 4.5;
% ISIs = 316/1000;
% AMPs = 38/10000;

Ca = 10.^(-pCas+6);

versions = {'parms_v4', 'parms_v1d'};
mcode = [1 1 1];

RMSD = nan(length(ISIs), length(AMPs), length(pCas), iFs(end));
[output_mainfolder, modelname, ~, ~] = get_folder_and_model(mcode);

for iF = iFs  
    for i = 1:length(Ca)
        for m = 1:length(AMPs)
            
            AMP = AMPs(m);
            
            for n = 1:length(ISIs)
                ISI = ISIs(n);
                
                for iii = 1:2
                    version = versions{iii};
                    
                    filename = [output_mainfolder{2}, '\', version, '\', modelname,'\',fibers{iF}, '\pCa=',num2str(pCas(i)*10), '\', fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'];
                    disp(filename)
                    
                    if exist(filename, 'file')
                        
                        try
                            load(filename, 'tis','Cas','vis','Lis','oFi','parms', 'ts')
                            data(iii).oFi = oFi;
                        catch
                            disp(['Unable to read: ', filename])
                        end
                        
                    else
                        if iii == 2
                            data(iii).oFi = nan(size(oFi));
                        else
                            keyboard
                        end
                    end
                        
                end
                
                RMSD(n,m,i,iF) = sqrt(mean((data(1).oFi-data(2).oFi).^2, 'omitnan'));                
            end
        end
    end
end
aAMPs = repmat(AMPs, length(ISIs),1);
aISIs = repmat(ISIs(:), 1, length(AMPs));

%%
disp(max(RMSD(:)))

close all
figure(1)

for i = 1:size(RMSD,3)
    nexttile
    surf(aAMPs, aISIs, RMSD(:,:,i,iF))
    zlim([0 max(RMSD(:))])
    set(gca, 'yscale', 'log')
end

% 
% %% find the worst trial
% k = 1;
% 
% [i,j] = find(RMSD(:,:,k,iF) == max(RMSD(:,:,k,iF),[],'all', 'omitnan'));
% 
% AMP = aAMPs(i,j)
% ISI = aISIs(i,j)
% pCa = pCas(k)

return
%% identify bad trials
clc

for ii = 1:size(RMSD,3)
    disp(['pCa = ', num2str(pCas(ii))])
    [i,j] = find(RMSD(:,:,ii,iF)>.03);

%     cd([output_mainfolder{2}, '\', version, '\', modelname,'\',fibers{iF}, '\pCa=',num2str(pCas(ii)*10)])

    for k = 1:length(i)

        AMP = aAMPs(i(k),j(k));
        ISI = aISIs(i(k),j(k));

        filename = [output_mainfolder{2}, '\', version, '\', modelname,'\',fibers{iF}, '\pCa=',num2str(pCas(ii)*10), '\', fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'];

        if exist(filename, 'file')
            delete(filename)
        end
        
        disp([filename, ' exist: ', num2str(exist(filename, 'file'))])
    end
    disp(' ')
end

%% 

