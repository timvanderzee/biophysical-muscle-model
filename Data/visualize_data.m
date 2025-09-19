clear all; close all; clc
[username, githubfolder] = get_paths();

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
cd(['C:\Users\',username,'\OneDrive - KU Leuven\9. Short-range stiffness\matlab\data'])
ISIs = [1 10 100 316 1000 3160 10000]/1000;

for i = 2:length(fibers)
    
    load([fibers{i},'_cor_new.mat'],'data')
    
    for n = 1:size(data.Fexp, 3)
        for m = 1:size(data.Fexp, 4)
 
                F = data.Fexp(:,:,n,m);
                L = data.Lexp(:,:,n,m);
                t = data.texp(:,:,n,m);
                
                figure(1)
                set(gcf,'name', [fibers{i}, '_ISI=',num2str(ISIs(n))])
                subplot(2,4,m)
                plot(t, F); 
                ylim([0 2])
                box off

                
%                 Fm(k,i) = mean(F(:1:1000));
            
%             figure(2)
%             subplot(3,4,i)
%             plot(data.pCas, Fm(:,i),'o')
%             yline(0.1,'k--')
        end
        drawnow
        pause
    end
end

