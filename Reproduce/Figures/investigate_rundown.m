
close all
for j = 1:11
        
    figure(j)
%     set(gcf, 'name', num2str(j))
        
    for i = 1:7
    
        nexttile
        histogram(fac(i,j,:,:));
        title(num2str(pCas(i)))
    end
end

for j = 1:11
    figure(j)
    
    nexttile

    plot(pCas, F_pCa.F0(j,:),'o'); hold on
    
    for n = 1:length(ISIs) % 3
        for m = 1:length(AMPs) % 7
            plot(pCas,  F_pCa.F02(:,j,n,m), '.'); hold on
        end
    end
    
    xlim([4 10])
end

%% determine average deviation
fac(1,[2 7],:,:) = nan;