clear all; close all; clc


load('SRS_data.mat', 'SRSrel', 'F0s', 'th')


%%
load('test_model_output_v6.mat','Stest', 'Scond', 'AMPs', 'iFs', 'pCas', 'ISIs', 'F0')
SRSrel_m = Stest./Scond(:,:,AMPs == .0383,:,2);
F0s_m = F0;

%%
ISIid = 1;
AMPid = 7;

close all
figure(1)
plot(F0s_m(:,ISIid, AMPid, 1, 2),SRSrel_m(:,ISIid, AMPid, 1, 2)); hold on
plot(F0s(:, AMPid,ISIid, 1), SRSrel(:,AMPid, ISIid, 1),'o')

iFs = [1 2 3, 5, 6, 7, 8, 10, 11];

for i = iFs
    SSE(i) = sum((SRSrel(:,AMPid, ISIid, i) - SRSrel_m(:,ISIid, AMPid, i, 2)).^2,'omitnan');
    SST(i) = sum((SRSrel(:,AMPid, ISIid, i) - mean(SRSrel(:,AMPid, ISIid, i),'omitnan')).^2,'omitnan');
end

%%
R2 = 1 - sum(SSE) ./ sum(SST)

%%

SST = sum((mean(SRSrel(:,AMPid, ISIid, :),4,'omitnan') - mean(mean(SRSrel(:,AMPid, ISIid, :),4,'omitnan'),'omitnan')).^2,'omitnan');
SSE = sum((mean(SRSrel(:,AMPid, ISIid, :),4,'omitnan') - mean(SRSrel_m(:,ISIid, AMPid, :,2),4,'omitnan')).^2,'omitnan');

R2 = 1 - sum(SSE) ./ sum(SST);

%% compute SSE from R2
R2s = [-24 .58];

SSE2 = -(R2s - 1) * sum(SST)

n = 4;

k = [5 10];
AIC = 2*k + n.*log(SSE2/n)