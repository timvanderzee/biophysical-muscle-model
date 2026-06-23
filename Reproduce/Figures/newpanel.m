clear all; close all; clc

cd('C:\Users\u0167448\Documents\GitHub\biophysical-muscle-model\Data')
load('force_pCa.mat')
close all
clc

%%
close all
pCas = [4.5 6.1:.1:6.4 6.6 9.0];

colors = parula(size(F0,1));

figure(1)
subplot(131)

for i = 1:size(F0,1)
plot(pCas, F0(i,:)./Fmax(i),'o','color', colors(i,:), 'markerfacecolor', colors(i,:), 'markersize', 3); hold on
end

legend('1',   '2',   '3',   '4',   '5',   '6',   '7',   '8',   '9',  '10',   '11')


for i = 1:length(pCas)
plot([pCas(i) pCas(i)], [0 1], 'k:'); hold on
end

for i = 1:size(F0,1)
plot(pCas, F0(i,:)./Fmax(i),'o','color', colors(i,:), 'markerfacecolor', colors(i,:), 'markersize', 3); hold on
end

legend('1',   '2',   '3',   '4',   '5',   '6',   '7',   '8',   '9',  '10',   '11')

xlabel('pCa')
ylabel('Activation level (F_0)')
box off


%%
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
k = 6;

iFs = 1:11;

T = .0843*2;

cd('C:\Users\u0167448\OneDrive - KU Leuven\9. Short-range stiffness\matlab\data')
load([fibers{iFs(k)},'_cor_new.mat'],'data')

t = data.texp - .003;
    
dts = [.01 .01];
cid = [1 3];

cs = [0 0 0; .5 .5 .5];

close all
for i = 1:2
    
    id = t(:,1,1,7)<.1 & t(:,1,1,7) > (-T-.1);
    id2 = t(:,1,1,7)<(-T+dts(i)) & t(:,1,1,7) > -T;

    id3 = t(:,1,1,7)<dts(i) & t(:,1,1,7) > 0;

    figure(1)
%     subplot(1,2,i)

    plot(data.Lexp(id,cid(i),1,7), data.Fexp(id,cid(i),1,7),':', 'marker', '.', 'color', cs(i,:), 'markersize', 8); hold on
    plot(data.Lexp(id2,cid(i),1,7), data.Fexp(id2,cid(i),1,7),'ro', 'markersize', 3)
    plot(data.Lexp(id3,cid(i),1,7), data.Fexp(id3,cid(i),1,7),'bo', 'markersize', 3)
    
    p1 = polyfit(data.Lexp(id2,cid(i),1,7), data.Fexp(id2,cid(i),1,7),1);
    p2 = polyfit(data.Lexp(id3,cid(i),1,7), data.Fexp(id3,cid(i),1,7),1);
    
    plot([0 .04], p1(2) + p1(1) * [0 .04], 'r--')
    plot([0 .04], p2(2) + p2(1) * [0 .04], 'b--')
    
    axis([0 .04 0 2])
    box off
    
    xlabel('Length')
    
    if i == 1
        ylabel('Force')
    end
    
%     if i == 2
%         subplot(1,3,3)
% 
%         plot(data.Lexp(id,cid(i),1,7), data.Fexp(id,cid(i),1,7),'k.'); hold on
%         plot(data.Lexp(id2,cid(i),1,7), data.Fexp(id2,cid(i),1,7),'ro','markerfacecolor', 'red', 'markersize', 3)
%         plot(data.Lexp(id3,cid(i),1,7), data.Fexp(id3,cid(i),1,7),'bo', 'markerfacecolor', 'blue', 'markersize', 3)
%     
%         plot([0 .04], p1(2) + p1(1) * [0 .04], 'r--')
%         plot([0 .04], p2(2) + p2(1) * [0 .04], 'b--')
% 
%         ylim([0    0.025])
%         box off
%         xlabel('Length')
%     end

end