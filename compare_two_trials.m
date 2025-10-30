clear all; close all; clc

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

iF = 2;

data1 = load(['C:\Users\u0167448\OneDrive - KU Leuven\9. Short-range stiffness\matlab\parms_v4\biophysical_full_alternative\', fibers{iF}, '\pCa=45\', fibers{iF},'_AMP=0_ISI=1.mat']);
data2 = load(['C:\Users\u0167448\OneDrive - KU Leuven\9. Short-range stiffness\matlab\parms_v2d\biophysical_full_alternative\', fibers{iF}, '\pCa=45\', fibers{iF},'_AMP=0_ISI=1.mat']);

plot(data1.tis, data1.oFi); hold on
plot(data2.tis, data2.oFi)
