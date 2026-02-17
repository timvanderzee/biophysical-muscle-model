clear all; close all; clc
[username, githubfolder] = get_paths();

iFs = [2,3,5,6,7,8,11];
    
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

AMPs = [0 12 38 121 216 288 383 682]/10000;
ISIs = [1 10 100 316 1000 3160 10000]/1000;
pCas = [4.5 6.1 6.2 6.3 6.4 6.6 9];

Ca = 10.^(-pCas+6);

versions = {'parms_v2d','parms_v4'};
mcode = [1 2 1];

RMSD = nan(length(ISIs), length(AMPs), length(pCas), iFs(end));
[output_mainfolder, modelname, ~, ~] = get_folder_and_model(mcode);

%%
for iF = iFs  
    for i = 1:length(Ca)
        for m = 1:length(AMPs)
            
            AMP = AMPs(m);
            
            for n = 1:length(ISIs)
                ISI = ISIs(n);
        
                filename1 = [output_mainfolder{2}, '\', versions{1}, '\', modelname,'\',fibers{iF}, '\pCa=',num2str(pCas(i)*10), '\', fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'];
                filename2 = [output_mainfolder{2}, '\', versions{2}, '\', modelname,'\',fibers{iF}, '\pCa=',num2str(pCas(i)*10), '\', fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'];

                if exist(filename1, 'file') && exist(filename2, 'file')
                    disp(filename1)

                    try
                        data1 = load(filename1, 'tis','oFi');
                        data2 = load(filename2, 'tis','oFi');
                        RMSD(n,m,i,iF) = sqrt(mean((data1.oFi(data1.tis>2)-data2.oFi(data1.tis>2)).^2, 'omitnan'));   
                         
                    catch
                        disp(['Unable to read: ', filename1])
                    end
                end            
            end
        end
    end
end

aAMPs = repmat(AMPs, length(ISIs),1);
aISIs = repmat(ISIs(:), 1, length(AMPs));

%%
disp(max(RMSD(:)))

close all

for iF = iFs
figure(iF)

for i = 1:size(RMSD,3)
    nexttile
    surf(aAMPs, aISIs, RMSD(:,:,i,iF)); hold on
    zlim([0 max(RMSD(:,:,:,iF), [], 'all')]); hold on
    
    surf(aAMPs, aISIs, .02 * ones(size(aAMPs)), 'facecolor', [.5 .5 .5])
    set(gca, 'yscale', 'log')
    
    % axis([0 .04 1e-3 1e1 0 max(RMSD(:,:,:,iF), [], 'all')])
end
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
% iF = 11;

for iF = iFs
    disp(fibers{iF})
    
for ii = 1:size(RMSD,3)
    disp(['pCa = ', num2str(pCas(ii))])
    [i,j] = find(RMSD(:,:,ii,iF)>.1);

%     cd([output_mainfolder{2}, '\', version, '\', modelname,'\',fibers{iF}, '\pCa=',num2str(pCas(ii)*10)])

    for k = 1:length(i)

        AMP = aAMPs(i(k),j(k));
        ISI = aISIs(i(k),j(k));

        filename = [output_mainfolder{2}, '\', version, '\', modelname,'\',fibers{iF}, '\pCa=',num2str(pCas(ii)*10), '\', fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'];

        if exist(filename, 'file')
%             delete(filename)
        end
        
        disp([filename, ' exist: ', num2str(exist(filename, 'file'))])
    end
    disp(' ')
end
end

%% 

