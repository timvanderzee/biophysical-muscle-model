clear all; close all; clc
cd(fileparts(which('reproduce_manuscript_figures.m')));
cd ..
cd ..
githubfolder = cd;

% get folders
cd('biophysical-muscle-model')
[datafolder, modelfolder] = get_paths(githubfolder);

figs = 5;

for fig = figs
    cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Reproduce', 'Figures'))

    %% Figures 2-4: Force traces    
    if fig > 1 && fig < 5
        Figs2to4(datafolder, modelfolder, githubfolder, fig)
    elseif fig == 5
        
        %% Figure 5: Force RMSD
        cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Reproduce', 'Figures'))
        Fig5(githubfolder)
        
    elseif fig == 6
        %% Figure 6: SRS
        cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Reproduce', 'Figures'))
        Fig6(githubfolder)
        
    elseif fig == 7
        %% Figure 7: SRS
        cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Reproduce', 'Figures'))
        Fig7(datafolder, modelfolder, githubfolder)
        
    elseif fig == 8
        %% Figure 8: Computational cost
        cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Reproduce', 'Figures'))
        Fig8(githubfolder)
    end
end

return
%% Export to PDF
cd('C:\Users\u0167448\OneDrive - KU Leuven\9. Short-range stiffness\revision\figures\PDF')

for i = figs
    figure(i)
    %     exportgraphics(gcf,['Fig', num2str(i), '.svg'],'ContentType','vector')
    exportgraphics(gcf,['FigS', num2str(i-1), '.pdf'],'ContentType','vector')
end