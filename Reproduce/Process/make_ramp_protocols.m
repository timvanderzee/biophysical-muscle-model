function[] = make_ramp_protocols(githubfolder, modelfolder, visualize)

save_results = 1;
redo = 0;

% conditions
AMPs = [0 12 38 121 216 288 383 532 682]/10000;
ISIs = [1 10 50 100 200 316 500 1000 3160 10000]/1000;
pCas = [4.5 6.1 6.2 6.3 6.4 6.6 9];
Ca = 10.^(-pCas+6);

if visualize
    if ishandle(1), close(1); end
    color = lines(1);
    
    figure(1)
    subplot(311)
    h(1) = line('xdata', [], 'ydata', [], 'color', color, 'linewidth', 1.5);
    xlim([0 14])
    ylim([0 35])
    box off
    title('Calcium')

    subplot(312)
    h(2) = line('xdata', [], 'ydata', [], 'color', color, 'linewidth', 1.5);
    xlim([0 14])
    ylim([0 .08])
    box off
    title('Fiber length')
    
    subplot(313)
    h(3) = line('xdata', [], 'ydata', [], 'color', color, 'linewidth', 1.5);
    xlim([0 14])
    ylim([0 .5])
    box off              
    title('Fiber velocity')
    xlabel('Time (s)')
end

for i = 1:length(Ca)

    output_foldername = fullfile(modelfolder, 'Protocols', ['pCa=',num2str(pCas(i)*10)]);
    
    if ~isfolder(output_foldername)
        mkdir(output_foldername)
    end
    
    for ii = 1:length(AMPs)
        
        AMP = AMPs(ii);
        dTt = .0383/.4545; % test stretch (= constant)
        dTc = AMP / .4545; % conditioning stretch
        
        for jj = 1:length(ISIs)
            ISI = ISIs(jj);
            
            filename = fullfile(output_foldername, ['AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat']);
            
            if ~exist(filename, 'file') || redo
                
                tiso = dTt*3+dTc*2+ISI + 2;
                
                dt = .001; % gives 10 points in SRS zone
                N = round(tiso / dt) + 1;
                
                cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Common', 'Functions'))
                [tis, Cas, Lis, vis, ts, Ts] = create_input(tiso, dTt, dTc, ISI, Ca(i), N);
                 
                if save_results
                    disp(['Saving ', filename]);
                    save(filename, 'tis','Cas','vis','Lis','Ts','ts')
                end
                
                if visualize
                    figure(1)
                    subplot(311)
                    set(h(1), 'xdata', tis, 'ydata', Cas)
   
                    
                    subplot(312)
                    set(h(2), 'xdata', tis, 'ydata', Lis)
    
                    
                    subplot(313)
                    set(h(3), 'xdata', tis, 'ydata', vis)

                    drawnow
                    
                end
            end
        end
    end
end



