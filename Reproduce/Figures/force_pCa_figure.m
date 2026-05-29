cd('C:\Users\u0167448\Documents\GitHub\biophysical-muscle-model\Data')
load('force_pCa.mat')
close all
clc

%%
close all
pCas = [4.5 6.1:.1:6.4 6.6 9.0];

colors = parula(size(F0,1));

figure(1)

for i = 1:length(pCas)
plot([pCas(i) pCas(i)], [0 1], 'k:'); hold on
end

for i = 1:size(F0,1)
plot(pCas, F0(i,:)./Fmax(i),'o','color', colors(i,:), 'markerfacecolor', colors(i,:), 'markersize', 3); hold on
end



xlabel('pCa')
ylabel('Activation level (F_0)')
box off

