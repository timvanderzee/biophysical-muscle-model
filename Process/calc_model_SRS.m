function[] = calc_model_SRS(githubfolder, modelfolder, iii, discretized_model)
visualize = 0;

% fibers
iFs = [2,3,5,6,7,8,11];
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

% models
mcodes = [2 2 1; 2 1 1; 1 1 3; 1 1 2; 1 1 1; 1 2 1];
mcode = mcodes(iii,:);

% conditions
AMPs = [0 12 38 121 216 288 383 532 682]/10000;
ISIs = [1 10 50 100 200 316 500 1000 3160 10000]/1000;
pCas = [4.5 6.1 6.2 6.3 6.4 6.6 9];
Ca = 10.^(-pCas+6);

%% folders
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
else
    modelname = [model,'_',mvar];
end

% folder containing the .mat files
fullmodelfolder = [modelfolder, subfolder];

% folder where output is saved
outputfolder = [githubfolder, '/biophysical-muscle-model/Model output/SRS/', subfolder];

%% loop over fibers, pCas, AMPs, ISIs
% pre-allocate
F0      = nan(length(pCas), length(ISIs),length(AMPs), 11);
Scond   = nan(length(pCas), length(ISIs),length(AMPs), 11);
Stest   = nan(length(pCas), length(ISIs),length(AMPs), 11);

for iF = iFs
    for i = 1:length(Ca)
        for ii = 1:length(AMPs)
            
            AMP = AMPs(ii);
            dTt = .0383/.4545; % test stretch (= constant)
            dTc = AMP / .4545; % conditioning stretch
            
            for jj = 1:length(ISIs)
                
                ISI = ISIs(jj);
                
                tiso = dTt*3+dTc*2+ISI + 2;

                filename = [fullmodelfolder, '\', modelname,'\',fibers{iF}, '\pCa=',num2str(pCas(i)*10), '\', fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'];
                disp(filename)
                
                if exist(filename, 'file')
                    
                    try
                        load(filename, 'tis','Cas','vis','Lis','oFi','parms', 'ts')
                        
                        cd([githubfolder, '\biophysical-muscle-model\Functions']) 
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
                    catch
                        keyboard
                        disp('Unable to read file')
                    end
                else
                    disp('File does not exist')
                end
            end
        end
    end
end


%% save
if ~isfolder(outputfolder)
    mkdir(outputfolder)
end

cd(outputfolder)
save([modelname, '_SRS.mat'])


