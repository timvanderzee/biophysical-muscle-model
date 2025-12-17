clear all; close all; clc

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

iF = 11;
AMP = 0;
ISI = 1;
pCa = 61;

model = 'biophysical_full_alternative';
% model = 'biophysical_full_regular';

data1 = load(['C:\Users\u0167448\OneDrive - KU Leuven\9. Short-range stiffness\matlab\parms_v6\',model, '\', fibers{iF}, '\pCa=', num2str(pCa), '\', fibers{iF},'_AMP=', num2str(AMP),'_ISI=', num2str(ISI), '.mat']);
data2 = load(['C:\Users\u0167448\OneDrive - KU Leuven\9. Short-range stiffness\matlab\parms_v4d\',model, '\', fibers{iF}, '\pCa=', num2str(pCa), '\', fibers{iF},'_AMP=', num2str(AMP), '_ISI=', num2str(ISI), '.mat']);



plot(data1.tis, data1.oFi); hold on
plot(data2.tis, data2.oFi)
