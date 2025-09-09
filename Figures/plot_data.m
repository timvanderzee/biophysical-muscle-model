
function[] = plot_data(texp, Lexp, Fexp, Tsrel)

for kk = 1:size(texp,2)
    tids = sort([Tsrel(kk,:) .15]);
    
    subplot(3,3,kk)
    plot(texp(:,kk), Lexp(:,kk)*100, 'color', [.5 .5 .5], 'linewidth', 2); hold on
    
    for ii = 1:(length(tids)-2)
        plot([tids(ii) tids(ii)], [0 200], 'k:'); hold on
    end
    
    axis([-.35 .25 0 5])
    box off
    
    if kk == 1
        ylabel('Length (%L_0)')
    end
    
    
    subplot(3,3,kk+3)
    plot(texp(:,kk), Fexp(:,kk)*100, 'color', [.5 .5 .5], 'linewidth', 2); hold on
    
    for ii = 1:(length(tids)-2)
        plot([tids(ii) tids(ii)], [0 200], 'k:'); hold on
    end
    
    axis([-.35 .25 0 max(Fexp(:,kk)*100)*1.2])
    box off
    
   
    
    if kk == 1
        ylabel('Force (%F_0)')
    end
        
end
end