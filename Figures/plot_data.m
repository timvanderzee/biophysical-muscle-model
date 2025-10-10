function[] = plot_data(texp, Lexp, Fexp, color, ls)

n = size(texp,2);

for kk = 1:n
    
    subplot(4,n,kk)
    plot(texp(:,kk), Lexp(:,kk)*100, 'color', color, 'linewidth', 2, 'linestyle',ls); hold on
    
    axis([-.35 .15 -.1 5])
    box off

    set(gca, 'Fontsize', 6)
    
    if kk == 1
       
        ylabel('Length (%L_0)', 'Fontsize', 8)
    end
    
    subplot(4,n,[kk+n kk+3*n])
    plot(texp(:,kk), Fexp(:,kk)*100, 'color', color, 'linewidth', 2, 'linestyle',ls); hold on
        
    axis([-.35 .15 0 max(Fexp(:,kk)*100)*1.2])
    box off
    
    set(gca, 'Fontsize', 6)
    
    if kk == 1
        
        ylabel('Force (%F_0)', 'Fontsize', 8)
    end
        
end
end