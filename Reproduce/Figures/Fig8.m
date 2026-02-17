function[] = Fig8(githubfolder)

cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Reproduce', 'Model output'))
load('output_effect_of_bins_on_cost.mat', 'tsim', 'Fi', 'bw', 'tvec', 'Ns', 'r')
titles = {'Fixed strain range (± 15 \epsilon_{ps})', ['Fixed bin width (', num2str(bw), ' \epsilon_{ps})']};

figure(8)

% color = get(gca,'colororder');
color = flip(copper(3));

Ntot = 1 * ones(size(tsim));

% Fint = mean(mean(Fi(:,end,1:2,:),4),3);
error = nan(size(Fi));

% Fint = Fi(:,end,1,1,1);
id = tvec > 1.9;
Fint = mean(mean(Fi(id,end,1,1:2,1), 4),5);

for f = 1:3
    error(id,:,:,f,:) = abs(Fi(id,:,:,f,:) - Fint(:));
end

eps = squeeze(mean(error, 'omitnan'));

for j = 1:2 % tests
    for f = 1:2 % methods
        
        subplot(2,3,(j-1)*3+1)
        set(gca, 'fontsize', 6)
        
        semilogx(Ns, median(tsim(:,j,f,:),4,'omitnan') ./ Ntot(:,j,f,r),'o--', 'color', color(f,:), ...
            'markerfacecolor', color(f,:), 'markersize', 4); hold on
        box off
        xlabel('Number of strain bins', 'fontsize', 8)
        ylabel('Time per sim. (s)', 'fontsize', 8)
        
        title(titles{j},'fontsize', 8)
        subtitle('Processing time', 'fontsize', 8)
        
        xlim([10 1.001e4])
        ylim([0 5])
        
        subplot(2,3,(j-1)*3+2)
        set(gca, 'fontsize', 6)
        semilogx(Ns, median(eps(:,j,f,:),4),'o--', 'color', color(f,:), ...
            'markerfacecolor', color(f,:), 'markersize', 4); hold on
        box off
        subtitle('Force error', 'fontsize', 8)
        xlabel('Number of strain bins', 'fontsize', 8)
        ylabel('Force error (F_0)', 'fontsize', 8)
        
        xlim([10 1.001e4])
        ylim([0 .1])
                        
        subplot(2,3,(j-1)*3+3)
        set(gca, 'fontsize', 6)
        plot(median(tsim(:,j,f,:),4,'omitnan') ./ Ntot(:,j,f), median(eps(:,j,f,:),4),'o--', 'color', color(f,:), ...
            'markerfacecolor', color(f,:), 'markersize', 4); hold on
        
        box off
        subtitle('Error versus time trade-off', 'fontsize', 8)
        xlabel('Time per simulation (s)', 'fontsize', 8)
        ylabel('Force error (F_0)', 'fontsize', 8)
        
        ylim([0 .1])
        xlim([0 2])
        %         xlim([Ns(1) Ns(end-1)])
        
    end
end

for i = [1,2,4,5]
    subplot(2,3,i)
%     yticks(10.^(-6:4));
    xticks(10.^(-6:4));
end

% plot approximated
f = 3;
figure(8)

for j = 1:2 % tests       
    subplot(2,3,(j-1)*3+1)
    yline(median(min(tsim(:,j,f,:),[],4),'all'),'--', 'color', color(f,:), 'linewidth', 1)

    subplot(2,3,(j-1)*3+2)
    yline(median(eps(:,j,f,:),'all'),'--', 'color', color(f,:), 'linewidth', 1)

    subplot(2,3,(j-1)*3+3)
    plot(median(min(tsim(:,j,f,:),[],4),'all'), median(eps(:,j,f,:),'all'),'ko', 'markerfacecolor', color(f,:), 'markersize', 4)
end

%%
figure(8)
% add second y-axis
N = length(tvec);

for j = 1:2
    subplot(2,3,(j-1)*3+1)
    
    yyaxis left
    
%     xlim([0 1e3])
    
    yl = get(gca, 'ylim');
    
    yyaxis right
    set(gca, 'yscale', 'log', 'ycolor', [.3 .3 .3])
    ylim(yl/N)
    ylabel('Time per eval. (s)','fontsize', 8)
    
    yticks(10.^(-10:3));

end



%%
figure(8)
set(gcf,'units','centimeters','position',[5 2 19 10])

subplot(231)
legend('Traditional method', 'Method of charac.', 'Gaussian approx.', 'location', 'best', 'fontsize', 6)


% subplot(326)
% legend('Tradit.', 'Charac.', 'Approx.', 'location', 'best')
% legend boxon

return
%% figure size
cd('C:\Users\u0167448\OneDrive\9. Short-range stiffness\figures\MAT')
figname = 'Fig10.png';
savefig = 1;
if savefig
    figure(8)
    exportgraphics(gcf,figname)
end
end
