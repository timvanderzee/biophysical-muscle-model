clear all; close all; clc


load('SRS_data.mat', 'SRSrel', 'F0s', 'th')


%%
load('test_model_output_v8.mat','Stest', 'Scond', 'AMPs', 'iFs', 'pCas', 'ISIs', 'F0')
% SRSrel_m = Stest./Scond(:,:,AMPs == .0383,:,2);
% F0s_m = F0;

% average 
iFs = [1 2 3, 5, 6, 7, 8, 10, 11];
th = [0 .05 .1 .25 .7 1.5];
SRSrel_m = nan(length(th)-1,7,8,length(iFs));
F0s_m = nan(length(th)-1,7,8,length(iFs));

for k = 1:length(iFs)
    for i = 1:length(th)-1
        id = F0(:,1,7,k) > th(i) & F0(:,1,7,k) <= th(i+1);
        
        SRSrel_m(i,:,:,k) = mean(Stest(id,:,:,k),1,'omitnan') ./ mean(Scond(id,1,7,k),'all', 'omitnan');
        F0s_m(i,:,:,k) = mean(F0(id,:,:,k), 1,'omitnan');
    end
end

%%
close all


for i = 1:3
    
    if i == 1
        ISIid = 1;
        pCaid = 1:size(SRSrel,1);
        AMPid = 7;
        x = pCas(pCaid);
    elseif i == 2
        ISIid = 1;
        pCaid = 4;
        AMPid = [1 2 3 4 7];
        x = AMPs(AMPid);
    elseif i == 3
        ISIid = [1, 3, 4, 5, 7];
        pCaid = 4;
        AMPid = 7;
        x = ISIs(ISIid);
    end
    
    SST(i) = sum((squeeze(mean(SRSrel(pCaid,ISIid, AMPid, :),4,'omitnan')) - mean(mean(SRSrel(pCaid,ISIid, AMPid, :),4,'omitnan'),'omitnan')).^2,'omitnan');
    SSE(i) = sum((squeeze(mean(SRSrel(pCaid,ISIid, AMPid, :),4,'omitnan')) - squeeze(mean(SRSrel_m(pCaid,ISIid, AMPid, :),4,'omitnan'))).^2,'omitnan');

    figure(1)
    subplot(1,3,i)
    plot(x, (squeeze(mean(SRSrel(pCaid,ISIid, AMPid, :),4,'omitnan')))); hold on
    plot(x, squeeze(mean(SRSrel_m(pCaid,ISIid, AMPid, :),4,'omitnan')))

end

R2 = 1 - SSE ./ SST;

return

%% compute SSE from R2
R2s = [-24 -.02 .58];
% R2s = [.08 .66 .72];
% R2s = [-1.2 .49 .73];

SSE2 = -(R2s - 1) * sum(SST)

n = 4;

k = [5 10 12];
AIC = 2*k + n.*log(SSE2/n)

% AICc = AIC + (2*k.^2 + 2*k) ./ (n-k-1)


