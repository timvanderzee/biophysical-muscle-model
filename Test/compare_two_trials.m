clear all; close all; clc

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

iF = 3;
AMP = 682;
ISI = 1;

data1 = load(['C:\Users\u0167448\OneDrive - KU Leuven\9. Short-range stiffness\matlab\parms_v4\biophysical_full_alternative\', fibers{iF}, '\pCa=45\', fibers{iF},'_AMP=', num2str(AMP),'_ISI=', num2str(ISI), '.mat']);
data2 = load(['C:\Users\u0167448\OneDrive - KU Leuven\9. Short-range stiffness\matlab\parms_v2d\biophysical_full_alternative\', fibers{iF}, '\pCa=45\', fibers{iF},'_AMP=', num2str(AMP), '_ISI=', num2str(ISI), '.mat']);

plot(data1.tis, data1.oFi); hold on
plot(data2.tis, data2.oFi)
