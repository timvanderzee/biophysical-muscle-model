clear all; close all; clc
[username, githubfolder] = get_paths();

Kss = 1:7;
tiso = 3;

iFs = 1:11;

F0 = nan(length(Kss), 7,8,length(iFs));
SRS_pre = nan(length(Kss), 7,8,length(iFs));
SRS_post = nan(length(Kss), 7,8,length(iFs));
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

visualize = 0;

for k = 1:length(iFs)
    cd(['C:\Users\',username,'\OneDrive - KU Leuven\9. Short-range stiffness\matlab\data'])
    load([fibers{iFs(k)},'_cor_new.mat'],'data')
    
    disp([fibers{iFs(k)}])
    
    for n = 1:7
        disp(n)
        for m = 1:8
            
            Data = prep_data(data,n,m,Kss,tiso);
            
            if ~isempty(Data.t)
                ts = 0:tiso:(tiso*(length(Kss)-1));
                
                % pre-allocate
                ns = nan(length(Kss), 2);
                os = nan(length(Kss), 2);
                ds = nan(length(Kss), 2);
                
                [id0,id1,id2] = get_indices(Data.t, tiso, ts, Data.dTt, Data.dTc, Data.ISI, Data.Ca(Kss));
                
%                 id0 = repmat(1:100,7,1);

                if visualize
                    figure(1)
%                     nexttile
                    plot(Data.t, Data.F,'.'); hold on


                    figure(2)
%                     nexttile
                    plot(Data.L, Data.F,'.'); hold on
                end
                
                for i = 1:length(Kss)
                    
                    % filter
%                     fs = 1000;
%                     fc = 100;
%                     Wn = fc / (.5*fs);
%                     [b,a] = butter(2, Wn);
%                     
%                     F = nan(size(Data.F));
%                     L = nan(size(Data.F));
%                     
%                     id = isfinite(Data.F);
%                     F(id) = filtfilt(b,a, Data.F(id));
%                     L(id) = filtfilt(b,a, Data.L(id));
                    
                    F = Data.F;
                    L = Data.L;
                    
                    dp1 = polyfit(L(id1(i,:)), F(id1(i,:)), 1);
                    
                    SRS_post(i,n,m,k) = dp1(1);

                    if visualize
                        figure(1)
                        plot(Data.t(id1(i,:)), F(id1(i,:)), 'r.')

                        figure(2)
                        plot(L(id1(i,:)), F(id1(i,:)), 'r.')
                    end

                    
                    % only for ISI = 1 ms
                    % only for AMP = 3.83%
                    if m == 7 && n < 5
                        dp2 = polyfit(L(id2(i,:)), F(id2(i,:)), 1);
                        
                        F0(i,n,m,k) = mean(F(id0(i,:)));
                        SRS_pre(i,n,m,k) = dp2(1);
   
                        if visualize
                            
                            figure(1)
                            plot(Data.t(id0(i,:)), F(id0(i,:)), 'g.')
                            plot(Data.t(id2(i,:)), F(id2(i,:)), 'y.')
                            hold off

                            figure(2)
                            plot(L(id0(i,:)), F(id0(i,:)), 'g.')
                            plot(L(id2(i,:)), F(id2(i,:)), 'y.')
                        end
                    end
                end
            end
        end
    end
end

return

%% average some pCas
th = [0 .05 .1 .25 .7 1.5];
SRSrel = nan(length(th)-1,7,8,length(iFs));
F0s = nan(length(th)-1, 7,8,length(iFs));

for k = 1:length(iFs)
    for i = 1:length(th)-1
        id = F0(:,1,7,k) > th(i) & F0(:,1,7,k) <= th(i+1);
        
        SRSrel(i,:,:,k) = mean(SRS_post(id,:,:,k),1,'omitnan') ./ mean(SRS_pre(id,:,7,k),'all', 'omitnan');
        F0s(i,:,:,k) = mean(F0(id,:,:,k), 1,'omitnan');
    end
end

%% save
cd([githubfolder, '\biophysical-muscle-model\Data'])
save('SRS_data_v3.mat', 'SRSrel', 'F0s', 'SRS_post', 'SRS_pre', 'th', 'F0')

