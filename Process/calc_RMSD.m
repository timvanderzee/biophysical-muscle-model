function[] = calc_RMSD(githubfolder, datafolder, modelfolder, iii, discretized_model)

% save folder
if iii > 2 % biophysical model
    if discretized_model 
        subfolder =  '\Discretized\';
    else
        subfolder = '\Approximated\';
    end
else % Hill-type model
    subfolder = '\Hill\';
end

save_output = 1;

Kss = 1:7;
tiso = 3;

iFs = [2,3,5,6,7,8,11];
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
mcodes = [2 2 1; 2 1 1; 1 1 3; 1 1 2; 1 1 1; 1 2 1];

visualize = 0;

version = 'parms';

%% calc RMSD
AMPs = [0 12 38 121 216 288 383 682]/10000;
ISIs = [1 10 100 316 1000 3160 10000]/1000;
pCas = [4.5 6.1 6.2 6.3 6.4 6.6 9];
Ca = 10.^(-pCas+6);


RMSD = nan(length(pCas), length(ISIs), length(AMPs), iFs(end), 7);

mcode = mcodes(iii,:);

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
    modelname = [model,'_',cvar,'_',mvar];
    
    if discretized_model
        fullmodelfolder = [modelfolder, '\Discretized\'];
    else
        fullmodelfolder = [modelfolder, '\Approximated\'];
    end
    
else
    modelname = [model,'_',mvar];
    fullmodelfolder = [modelfolder, '\Hill\'];
end

for iF = iFs
    
    years = {'2017', '2018'};
    for m = 1:length(years)
        if contains(fibers{iF}, years{m})
            fullfolder = [datafolder, '\', years{m}];
        end
    end
    
    cd(fullfolder)
    load([fibers{iF},'.mat'],'data')
    
    for i = 1:length(Ca)
        cd(fullmodelfolder)
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
                        
                        mD = AMP == data.AMPs/10000;
                        nD = ISI == data.ISIs/1000;
                        iD = round(pCas(i),1) == round(data.pCas,1);
                        
                        texp = data.texp(:,iD,nD,mD) - .005;
                        Fexp = data.Fexp(:,iD,nD,mD);
                        
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
                                RMSDs = sqrt((oFii(id) - Fexp(id)).^2) * 100;
                                RMSD(i,n,m,iF,ii) = sqrt(mean((oFii(id) - Fexp(id)).^2, 'omitnan'));
                                
                                if visualize
                                    disp(RMSD(i,n,m,iF,ii))
                                    figure(1)
                                    plot(tm, oFi, '-', texp, Fexp, '-', texp(id), oFii(id), '.')
                                    pause
                                end
                                
                                
                            else
                                RMSD(i,n,m,iF,ii) = nan;
                                %                                     disp('Does not exist')
                            end
                            
                            
                            if RMSD(i,n,m,iF,ii) > .2
                                
                                if i == 1 && iF == 2 && (AMP < .03 || ISI > 1)
                                    disp('bad trial')
                                    
                                else
                                    %                                         figure(iF)
                                    %                                         nexttile
                                    %                                         plot(tm, oFi, '-', texp, Fexp, '-', texp(id), oFii(id), '.')
                                    %                                         title([filename, ' - pCa=',num2str(pCas(i)*10)])
                                end
                                %                                     keyboard
                            end
                        end
                        
                    catch
                        disp(['Unable to read: ', filename])
                    end
                    
                else
                    disp('Does not exist')
                end
                
            end
        end
    end
end


%% save

if save_output
    fullmodelfolder = [githubfolder, '/biophysical-muscle-model/Model output/RMSD/', subfolder];
    
    if ~isfolder(fullmodelfolder)
        mkdir(fullmodelfolder)
    end
    
    cd(fullmodelfolder)
    
    save([modelname, '_RMSD.mat'])
end


