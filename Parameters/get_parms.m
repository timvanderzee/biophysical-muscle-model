clear all; close all; clc
% gets parameters from OneDrive and saves them
[username, githubfolder] = get_paths();

iFs = [1 2 3, 5, 6, 7, 8, 10, 11];

mcodes = [1 1 1; 1 1 3; 1 2 1; 2 1 1];

for i = 1:size(mcodes,1)
    % pparms = struct();
    clear pparms
    mcode = mcodes(i,:);
    
    for iF = iFs
        
        if mcode(2) == 2
            vs = '\full\';
        else
            vs = '\';
        end
        
        cd([githubfolder, '\muscle-thixotropy\new_model\get_variable'])
        [output_mainfolder, filename, opt_type, ~] = get_folder_and_model(mcode);
        
        disp(filename)
        output_folder = [opt_type,'\normalized\with_PE_optimized\2_trials'];
        
        output_dir = [output_mainfolder{1}, '\', filename,vs, output_folder];
        cd(output_dir)
        
        % get other parameters from other fiber
        % ii = 6; % fiber from which parameters are obtained
        load([filename,'_F', num2str(iF),'_best.mat'],'parms','exitflag','fopt','C0','Cbounds','model','P0','P')
        C = p_to_c(P, Cbounds);
        parms = C_to_parms(C, parms, parms.optvars);
        parms = calc_dependent_parms(parms);
        
        parms.act = 1;
        parms.Noverlap = 1;
        parms.vF_func =  @(vcerel,parms)parms.e(1)*log((parms.e(2)*vcerel./parms.vmax+parms.e(3))+sqrt((parms.e(2)*vcerel./parms.vmax+parms.e(3)).^2+1))+parms.e(4);

        pparms(iF) = parms;
    end
    
    cd([githubfolder, '\biophysical-muscle-model\Parameters'])
    save(['parms_',filename,'.mat'], 'pparms')
end