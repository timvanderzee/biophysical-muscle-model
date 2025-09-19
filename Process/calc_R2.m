clear all; close all; clc

load('SRS_data.mat', 'F0', 'SRS_post', 'SRS_pre')
F0d = F0;
%%
filenames = {'biophysical_full_regular_SRS.mat', 'Hill_regular_SRS.mat'};

for ff = 1:length(filenames)
    load(filenames{ff},'Stest', 'Scond', 'AMPs', 'pCas', 'ISIs', 'F0')
    % load('test_model_output_v8.mat','Stest', 'Scond', 'AMPs', 'iFs', 'pCas', 'ISIs', 'F0')
    % SRSrel_m = Stest./Scond(:,:,AMPs == .0383,:,2);
    % F0s_m = F0;

    % average 
%     iFs = [1 2 3, 5, 6, 7, 8, 10, 11];
    iFs = 1:11;
%     iFs = [2,3,5,6,7,8,11];
    th = [0 .07 .25 .7 1.5];
    SRSrel_m = nan(length(th)-1,7,8,iFs(end));
    F0s_m = nan(length(th)-1,7,8,iFs(end));
    
    SRSrel = nan(length(th)-1,7,8,iFs(end));
    F0s = nan(length(th)-1,7,8,iFs(end));
    
    for k = iFs
        for i = 1:length(th)-1
            id = F0(:,1,7,k) > th(i) & F0(:,1,7,k) <= th(i+1);

            SRSrel_m(i,:,:,k) = mean(Stest(id,:,:,k),1,'omitnan') ./ mean(Scond(id,1,7,k),'all', 'omitnan');
            F0s_m(i,:,:,k) = mean(F0(id,:,:,k), 1,'omitnan');
            
            id = F0d(:,1,7,k) > th(i) & F0d(:,1,7,k) <= th(i+1);        
            SRSrel(i,:,:,k) = mean(SRS_post(id,:,:,k),1,'omitnan') ./ mean(SRS_pre(id,1,7,k),'all', 'omitnan');
            F0s(i,:,:,k) = mean(F0d(id,:,:,k), 1,'omitnan');
        end
    end

    %% Compute R2 for specific conditions
%     close all
    for i = 1:3

        if i == 1
            ISIid = 1;
            pCaid = 1:size(SRSrel,1);
            AMPid = 7;
            x = mean(F0s(pCaid, ISIid, AMPid,:),4,'omitnan');
        elseif i == 2
            ISIid = 1;
            pCaid = 3;
            AMPid = [1  3 4 7];
            x = AMPs(AMPid);
        elseif i == 3
            ISIid = [1, 3, 4, 5, 7];
            pCaid = 3;
            AMPid = 7;
            x = ISIs(ISIid);
        end

        SST(i) = sum((squeeze(mean(SRSrel(pCaid,ISIid, AMPid, :),4,'omitnan')) - mean(mean(SRSrel(pCaid,ISIid, AMPid, :),4,'omitnan'),'omitnan')).^2,'omitnan');
        SSE(ff,i) = sum((squeeze(mean(SRSrel(pCaid,ISIid, AMPid, :),4,'omitnan')) - squeeze(mean(SRSrel_m(pCaid,ISIid, AMPid, :),4,'omitnan'))).^2,'omitnan');

        figure(ff)
        subplot(1,3,i)
        plot(x, (squeeze(mean(SRSrel(pCaid,ISIid, AMPid, :),4,'omitnan'))),'o'); hold on
        plot(x, squeeze(mean(SRSrel_m(pCaid,ISIid, AMPid, :),4,'omitnan')),'o')
        
        if i == 3
            set(gca,'XScale', 'log')
        end

    end
end

   R2 = 1 - SSE ./ SST;
return
%% compute R2 for all conditions
R2 = nan(size(SRSrel,1),1);

for i = 1:size(SRSrel,1)
SST = sum((mean(SRSrel(i,:,:,:),4,'omitnan') - mean(SRSrel(i,:,:,:),'all','omitnan')).^2,'all','omitnan');
SSE = sum((SRSrel(i,:,:,iFs) - SRSrel_m(i,:,:,:)).^2,'all','omitnan');

R2(i) = 1 - SSE ./ SST;
end

%%
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

return

%% compute SSE from R2
% R2s = [-24 -.02 .58];
% R2s = [.08 .66 .72];
% R2s = [-1.2 .49 .73];

% SSE2 = -(R2 - 1) .* SST

n = 5;
k = [10 10 10; 5 5 5];
AIC = 2*k + n.*log(SSE/n)

% AICc = AIC + (2*k.^2 + 2*k) ./ (n-k-1)


