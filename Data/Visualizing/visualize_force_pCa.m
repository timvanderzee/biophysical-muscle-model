clear all; close all; clc

load('active_trials', 'Fm')

pCas = [4.5 6.1 6.2 6.3 6.4 6.6 9];

for i = 1:size(Fm,2)
    nexttile
    plot(pCas, Fm(:,i), 'o-')
    yline(.1, 'k--')
    yline(.25, 'k--')
    yline(.7, 'k--')
end