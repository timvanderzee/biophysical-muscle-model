function[] = reorganize_data(datafolder, outputfolder)

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

pCas = [4.5 6.1:.1:6.4 6.6 9];
AMPs = [0 12 38 121 216 288 383 682];
ISIs = [1 10 100 316 1000 3162 10000];

for k = 1:length(fibers)
    
    cd(outputfolder)
    if ~exist([fibers{k},'_reorganized.mat'], 'file') % only if the file does not exist yet
        
        % initialize
        Force   = nan(2e4, length(pCas), length(ISIs), length(AMPs));
        Length  = nan(2e4, length(pCas), length(ISIs), length(AMPs));
        Time    = nan(2e4, length(pCas), length(ISIs), length(AMPs));
        
        pCa     = nan(length(pCas), length(ISIs), length(AMPs));
        cAMP    = nan(length(pCas), length(ISIs), length(AMPs));
        tAMP    = nan(length(pCas), length(ISIs), length(AMPs));
        ISI     = nan(length(pCas), length(ISIs), length(AMPs));
        F0      = nan(length(pCas), length(ISIs), length(AMPs));
        
        for i = 1:length(pCas)
            cd(datafolder)
            filename = [fibers{k},'_pCa_',strrep(num2str(pCas(i),'%.1f'),'.','pt'),'.mat'];
            disp(filename)
            
            if exist(filename, 'file')
                load(filename)
                
                for n = 1:length(ISIs)
                    for m = 1:length(AMPs)
                        
                        % search for data
                        f = 0;
                        for j = 1:length(new_data)
                            file_info = new_data(j).file_info_string;
                            fgain = str2double(file_info(20:21));
                            
                            
                            if new_data(j).testAmp == 383 && new_data(j).condVel == 455 && new_data(j).testVel == 455 && new_data(j).ISI == ISIs(n) && new_data(j).condAmp == AMPs(m)
                                
                                N = length(new_data(j).force);
                                Force(1:N,i,n,m) = new_data(j).force/1000 / fgain;
                                Length(1:N,i,n,m) = new_data(j).fibre_length;
                                Time(1:N,i,n,m) = new_data(j).time - new_data(j).Test_Onset/1000;
                                
                                pCa(i,n,m) = new_data(j).pCa;
                                cAMP(i,n,m) = new_data(j).condAmp;
                                tAMP(i,n,m) = new_data(j).testAmp;
                                ISI(i,n,m) = new_data(j).ISI;
                                id = Time(1:N,i,n,m) < (-ISI(i,n,m) - 2*(383/455));
                                F0(i,n,m) = mean(Force(id,i,n,m));
                                
                                f = 1;
                                %                 disp(length(Time))
                                break
                            else
                                
                                %             disp('Nothing found')
                            end
                        end
                        
                        if f
                            disp([fibers{k}, ' - pCa = ',num2str(pCas(i)), ' - ISI = ',num2str(ISIs(n)), ' - AMP = ',num2str(AMPs(m)), ': Found'])
                        else
                            disp([fibers{k}, ' - pCa = ',num2str(pCas(i)), ' - ISI = ',num2str(ISIs(n)), ' - AMP = ',num2str(AMPs(m)), ': Not found'])
                        end
                    end
                end
            end
        end

        %% save
        cd(outputfolder)
        save([fibers{k},'_reorganized.mat']);
    end
end
end