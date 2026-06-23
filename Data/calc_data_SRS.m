function[] = calc_data_SRS(current_folder, normalized_datafolder)

Kss = 1:7;
tiso = 3;

iFs = 1:11;

F0 = nan(length(Kss), 7,8,length(iFs));
SRS_pre = nan(length(Kss), 7,8,length(iFs));
SRS_post = nan(length(Kss), 7,8,length(iFs));
R2_pre = nan(length(Kss), 7,8,length(iFs));
R2_post = nan(length(Kss), 7,8,length(iFs));

LSRS_post = nan(20, 7,8,length(Kss),length(iFs));
FSRS_post = nan(20, 7,8,length(Kss),length(iFs));
LSRS_pre = nan(20, 7,8,length(Kss),length(iFs));
FSRS_pre = nan(20, 7,8,length(Kss),length(iFs));

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

visualize = 0;

tshifts = .001:.001:.005;

for tt = 3 %1:length(tshifts)

for k = iFs %1:length(iFs)
    cd(normalized_datafolder)
    load([fibers{iFs(k)},'_cor_new.mat'],'data')
    
    disp([fibers{iFs(k)}])
    
    for n = 1:7
        disp(n)
        for m = 1:8
            
            Data = prep_data_v2(data,n,m,Kss,tiso, tshifts(tt));
            
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
                    plot(Data.t, Data.F,'k.'); hold on


                    figure(2)
%                     nexttile
                    plot(Data.L, Data.F,'k.'); hold on
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
                    
                    LSRS_post(1:length(id1),n,m,i,k) = L(id1(i,:));
                    FSRS_post(1:length(id1),n,m,i,k) = F(id1(i,:));
                    LSRS_pre(1:length(id1),n,m,i,k) = L(id2(i,:));
                    FSRS_pre(1:length(id1),n,m,i,k) = F(id2(i,:));

                    dp1 = polyfit(L(id1(i,:)), F(id1(i,:)), 1);
                    
                    pred1 = polyval(dp1, L(id1(i,:)));
                    R2_post(i,n,m,k) = 1 - sum((pred1-F(id1(i,:))).^2) ./ sum((pred1-mean(pred1)).^2);
                    SRS_post(i,n,m,k) = dp1(1);
                    
                    
                    if R2_post(i,n,m,k) < 0.5 && mean(F(id1(i,:))) > 0.15
                        figure(10)
                        nexttile
                        plot(Data.t, F,'.')
%                         title(['
                    end
 

                    if visualize
                        figure(1)
                        plot(Data.t(id1(i,:)), F(id1(i,:)), 'r.')

                        figure(2)
                        plot(L(id1(i,:)), F(id1(i,:)), 'r.')
                        pause
                    end

                    
                    % only for ISI = 1 ms
                    % only for AMP = 3.83%
                    if m == 7 && n < 5
                        

                    
                        dp2 = polyfit(L(id2(i,:)), F(id2(i,:)), 1);
                        
                        pred2 = polyval(dp2, L(id2(i,:)));
                        
                        R2_pre(i,n,m,k) = 1 - sum((pred2-F(id2(i,:))).^2) ./ sum((pred2-mean(pred2)).^2);
                        
                        F0(i,n,m,k) = mean(F(id0(i,:)));
                        SRS_pre(i,n,m,k) = dp2(1);
   
                        if visualize
                            
                            figure(1)
                            plot(Data.t(id0(i,:)), F(id0(i,:)), 'g.')
                            plot(Data.t(id2(i,:)), F(id2(i,:)), 'm.')
                            hold off

                            figure(2)
                            plot(L(id0(i,:)), F(id0(i,:)), 'g.')
                            plot(L(id2(i,:)), F(id2(i,:)), 'm.')
                            
                            pause
                        end
                    end
                end
            end
        end
    end
end



%% average some pCas
% th = [0 .05 .1 .25 .7 1.5];
th = [0 .05 .3 .7 1.2];

SRSrel = nan(length(th)-1,7,8,length(iFs));
F0s = nan(length(th)-1, 7,8,length(iFs));

R2s_pre = nan(length(th)-1, 7,8,length(iFs));
R2s_post = nan(length(th)-1, 7,8,length(iFs));

SRSrel2 = nan(length(th)-1,7,8,length(iFs));
R2_2 = nan(length(th)-1,7,8,length(iFs));

for k = 1:length(iFs)
    for i = 1:length(th)-1
        id = F0(:,1,7,k) > th(i) & F0(:,1,7,k) <= th(i+1);
        
        SRSrel(i,:,:,k) = mean(SRS_post(id,:,:,k),1,'omitnan') ./ mean(SRS_pre(id,:,7,k),'all', 'omitnan');
        F0s(i,:,:,k) = mean(F0(id,:,:,k), 1,'omitnan');
        
        R2s_pre(i,:,:,k) = mean(R2_pre(id,:,:,k), 1, 'omitnan');
        R2s_post(i,:,:,k) = mean(R2_post(id,:,:,k), 1, 'omitnan');

        LSRSs_ref = mean(mean(LSRS_pre(:,1:4,7,id,k),2,'omitnan'), 4, 'omitnan');
        FSRSs_ref =  mean(mean(FSRS_pre(:,1:4,7,id,k),2,'omitnan'), 4, 'omitnan');
                    
        pref = polyfit(LSRSs_ref, FSRSs_ref,1);
                    
        for n = 1:7
            for m = 1:8
                LSRSs = mean(LSRS_post(:,n,m,id,k), 4, 'omitnan');
                FSRSs = mean(FSRS_post(:,n,m,id,k), 4, 'omitnan');

                if sum(isfinite(LSRSs)) > 0

                    ptest = polyfit(LSRSs, FSRSs,1);
                    
                    SRSrel2(i,n,m,k) = ptest(1)/pref(1);
                    
                    pred2 = polyval(ptest, LSRSs);

                    R2_2(i,n,m,k) = 1 - sum((pred2-FSRSs).^2) ./ sum((pred2-mean(pred2)).^2);
                        
                end
            end
        end
        

        
    end
end

%% remove bad data
R2_2(end,2:end,1:6,2) = nan;
R2_2(end,end,7,2) = nan;

SRSrel2(end,2:end,1:6,2) = nan;
SRSrel2(end,end,7,2) = nan;

%% save
cd(current_folder)
save(['SRS_data_DT20.mat'], 'SRSrel', 'F0s', 'SRS_post', 'SRS_pre', 'th', 'F0',  'SRSrel2')

end
end