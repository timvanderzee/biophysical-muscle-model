clear all; close all; clc


fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
cd('C:\Users\timvd\OneDrive - KU Leuven\9. Short-range stiffness\matlab\data')

for i = 1:length(fibers)

    load([fibers{i},'_cor_new.mat'],'data')

    n = 1;
    m = 1;
    Ks = 1:7;

    Data.F = [];
    Data.L = [];
    Data.t = [];
    Data.Ca = [];

    dTt = .0383/.4545; % test stretch (= constant)

    for k = 1:length(Ks)

        F = data.Fexp(1:1000,Ks(k),n,m);
        L = data.Lexp(1:1000,Ks(k),n,m);
        t = data.texp(1:1000,Ks(k),n,m);
        
        figure(1)
        subplot(3,4,i)
        plot(t, F); hold on
        ylim([0 2])
        box off
        
        Fm(k,i) = mean(F);
        
    end
    
    figure(2)
    subplot(3,4,i)
    plot(data.pCas, Fm(:,i),'o')
    yline(0.1,'k--')
end


