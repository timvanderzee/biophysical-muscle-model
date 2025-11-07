clear all; close all; clc

% paths
[username, githubfolder] = get_paths();

% binary inputs
save_results = 0;
redo = 0;
visualize = 0;
discretized_model = 0;

% versions
output_version = '_v5';
parms_version = '_v3'; 

% model
mcodes = [1 2 1];

% fibers
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
iFs = [2,3,5,6,7,8,11];

% conditions
AMPs = [0    0.0012    0.0038    0.0121    0.0216    0.0288    0.0383    0.0532    0.0682];
ISIs = [ 0.0010    0.0100    0.0500    0.1000    0.2000    0.3160    0.5000    1.0000    3.1600   10.0000];
pCas = [4.5 6.1 6.2 6.3 6.4 6.6 9];
Ca = 10.^(-pCas+6);

for iii = 1:size(mcodes,1)
    
    % load parameters
    mcode = mcodes(iii,:);
    [output_mainfolder, modelname, ~, ~] = get_folder_and_model(mcode);
    
    for iF = iFs
        
        % disp(filename)
        input_foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iF}];
        cd(input_foldername)
        load(['parms_',modelname, parms_version, '.mat'], 'newparms')
        allparms(iF) = update_parms(newparms);
    end
    
    parfor k = 1:length(iFs)
        iF = iFs(k);
        parms = allparms(iF);
        
        if contains(modelname, 'Hill')
            x0 = 0;
        elseif discretized_model
            
            n0 = zeros(size(parms.xi));
            x0 = [n0'; parms.x0(4:end)'];
            
        else
            x0 = parms.x0';
        end

        xp0 = zeros(size(x0));        
        
        for i = 1:length(Ca)
            X0 = x0;
            
            output_foldername = [output_mainfolder{2}, '\parms', output_version, '\', modelname,'\',fibers{iF}, '\pCa=',num2str(pCas(i)*10)];
            
            if ~isfolder(output_foldername)
                mkdir(output_foldername)
            end
               
            for ii = 1:length(AMPs)
                
                AMP = AMPs(ii);
                dTt = .0383/.4545; % test stretch (= constant)
                dTc = AMP / .4545; % conditioning stretch
                
                % when there is pre-stretch, make sure we don't miss it
                dtmax = max([dTc/2, 5e-3]);
                
                % if there is no pre-stretch, we can't miss it either
                if AMP == 0
                    dtmax = dTt;
                end
                
                dtmax = 1e-2;
                odeopt = odeset('maxstep',dtmax);
                
                for jj = 1:length(ISIs)
                    ISI = ISIs(jj);
                    
                    filename = [output_foldername, '\', fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'];
                    
                    if ~exist(filename, 'file') || ~save_results || redo
                        disp(filename);
                        tiso = dTt*3+dTc*2+ISI + 2;
                        
                        dt = .001; % gives 10 points in SRS zone
                        N = round(tiso / dt);
                        
                        [tis, Cas, Lis, vis, ts, Ts] = create_input(tiso, dTt, dTc, ISI, Ca(i), N);
                        
                        parms.ti = tis;
                        parms.vts = vis;
                        parms.Cas = mean(Cas);
                        parms.Lts = Lis * parms.gamma;
                        
                        % run simulation
                        if contains(modelname, 'Hill')
                            % simulate
                            sol = ode15i(@(t,y,yp) hill_type_implicit_v2(t,y,yp, parms), [0 max(tis)], X0, xp0, odeopt);
                            
                            % get SE length
                            Lse = parms.Lts - interp1(sol.x, sol.y(1,:), tis);
                            
                            % get force
                            if pCas(i) == 9
                                Fact = 0;
                            else
                                Fact = parms.Fse_func(Lse, parms) * parms.Fscale;
                            end
                            
                            oFi = Fact + parms.Fpe_func(parms.Lts, parms);
                            
                        else
                            
                            aTs = [0; Ts];
                            X0 = x0;
                            XP0 = xp0;
                            
                            vts = [0 .4545 -.4545 0 .4545 0 0];
                            
                            % splitting it up makes things much faster
                            tall = [];
                            Fall = [];
                            F2all = [];
                            Lall = [];
                            
                            % interval needs to have finite duration
                            nzi = find(diff(aTs) > 0);
                            
%                             odeopt = odeset('maxstep', 1e-3);
                            odeopt = [];
                            
                            for p = 1:(length(nzi)-1)
%                                 disp(p)
                                
                                % simulate
                                if discretized_model
                                    sol = ode15s(@(t,y,yp) fiber_dynamics_explicit_no_tendon_full(t,y, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, odeopt);
                                else
                                    sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, XP0, []);
                                end
                                
                                [~,xdot] = deval(sol, sol.x);
                                
                                X0 = sol.y(:,end);
                                XP0 = xdot(:,end);
                                
                                % get force
                                t = sol.x;
                                
                                if discretized_model
                                    L = sol.y(end-3,:);
                                    dlse = interp1(parms.ti, parms.Lts, t) - L;
%                                     F = parms.Fse_func(dlse, parms);
%                                     
                                    F = nan(1, length(sol.x));
                                    for iiii = 1:length(sol.x)
                                        n = sol.y(1:end-4,iiii);
                                        xi = parms.xi + L(iiii);
                                        F(iiii) = trapz(xi, xi .* n') + trapz(xi, n');
                                    end
                                    
                                else
                                    F = (sol.y(1,:) + sol.y(2,:));
                                end
                                
                                tall = [tall t];
                                Fall = [Fall F];
%                                 F2all = [F2all F2];
                                
                                
                                %                         figure(100)
                                %                         plot(t, F); hold on
                            end
                            
                            % find unique values
                            [~, ui] = unique(tall);
                            
                            % interpolate force
                            oFi = interp1(tall(ui), Fall(ui), tis) * parms.Fscale + parms.Fpe_func(parms.Lts, parms);
                        end
                        %                 toc
                        %%
%                         figure(1)
%                         plot(tall, Fall)
                        %%
                        % steady state
                        xs = sol.y(:,end);
                        
                        if save_results
                            cd(output_foldername);
%                             save(filename, 'tis','Cas','vis','Lis','oFi','parms','ts')
                            save_model_forces(filename, tis,Cas,vis,Lis,oFi,parms,ts)
                        end
                        
                        if visualize
                            figure(1)
                            subplot(411)
                            plot(tis, Cas, 'b', 'linewidth',1);
                            xlim([0 14])
                            
                            subplot(412)
                            plot(tis, vis,'b',  'linewidth',1);
                            xlim([0 14])
                            
                            subplot(413)
                            plot(tis, Lis,'b',  'linewidth',1);
                            xlim([0 14])
                            
                            subplot(414);
                            plot(tis, oFi,'b')
                            xlim([0 14])
                            
                            drawnow
                        end
                    end
                end
            end
        end
    end
end

