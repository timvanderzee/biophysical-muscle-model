function[] = calc_model_forces(githubfolder, modelfolder, iii, iFs, discretized_model, visualize)
save_results = 1;
redo = 1;
% visualize = 0;

% fibers
% iFs = [2,3,5,6,7,8,11];
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

% models
mcodes = [2 2 1; 2 1 1; 1 1 3; 1 1 2; 1 1 1; 1 2 1];
mcode = mcodes(iii,:);

% conditions
AMPs = [0 12 38 121 216 288 383 532 682]/10000;
ISIs = [1 10 50 100 200 316 500 1000 3160 10000]/1000;
pCas = [4.5 6.1 6.2 6.3 6.4 6.6 9];
Ca = 10.^(-pCas+6);

% folder containing the .mat files
[savefolder, ~, modelname] = get_model_folder(modelfolder, mcode, discretized_model);

%% loop over fibers, pCas, AMPs, ISIs
for iF = iFs
    
    % disp(filename)
    input_foldername = fullfile(githubfolder, 'biophysical-muscle-model', 'Reproduce', 'Parameters',fibers{iF});
    cd(input_foldername)
    load(['parms_',modelname, '.mat'], 'newparms')
    
    cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Reproduce', 'Parameters'))
    allparms(iF) = update_parms(newparms);
end

for k = 1:length(iFs)
    iF = iFs(k);
    parms = allparms(iF);
    
    if contains(modelname, 'Hill')
        x0 = 0;
        
    elseif discretized_model
        
        n0 = zeros(size(parms.xi));
        
        % find L0
        costfunc = @(L, Fce, parms) (Fce + (parms.kpe*(-L-parms.Lce0)) - (parms.kse0*(exp(parms.kse*L)-1))) .^2;
        dlse = fminsearch(@(L) costfunc(L, 0, parms), 0);
        lce0 = -dlse;
        
        x0 = [n0'; parms.x0(4:end)'];
        x0(end-3) = lce0;
        
        if isequal(mcode, [1 1 2]) || isequal(mcode, [1 1 3])
            x0(end-1) = 1;
        end
        
    else
        
        x0 = parms.x0';
        
        
        x0 = zeros(1,6);
        Q0 = .2;
        p0 = -1;
        q0 = .1;
        [Q00, Q20, lce0, Q10, Fse0] = find_steady_state(Q0, p0, q0, newparms, 'regular');
        x0 = [Q00 Q20 Fse0 0 0 0];
        
        if isequal(mcode, [1 1 2]) || isequal(mcode, [1 1 3])
            x0(5) = 1-Q00;
        end
    end
    
    xp0 = zeros(size(x0));
    
    for i = 1:length(Ca)
        X0 = x0;
        
        output_foldername = fullfile(savefolder, modelname,fibers{iF}, ['pCa=',num2str(pCas(i)*10)]);
        
        if ~isfolder(output_foldername)
            mkdir(output_foldername)
        end
        
        for ii = 1:length(AMPs)
            
            AMP = AMPs(ii);
            dTt = .0383/.4545; % test stretch (= constant)
            dTc = AMP / .4545; % conditioning stretch
            
            dtmax = 1e-2;
            odeopt = odeset('maxstep',dtmax);
            
            for jj = 1:length(ISIs)
                ISI = ISIs(jj);
                
                filename = fullfile(output_foldername, [fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat']);
                
                if ~exist(filename, 'file') || ~save_results || redo
                   
                    tiso = dTt*3+dTc*2+ISI + 2;
                    
                    dt = .001; % gives 10 points in SRS zone
                    N = round(tiso / dt);
                    
                    cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Common', 'Functions'))
                    [tis, Cas, Lis, vis, ts, Ts] = create_input(tiso, dTt, dTc, ISI, Ca(i), N);
                    
                    parms.ti = tis;
                    parms.vts = vis;
                    parms.Cas = mean(Cas);
                    parms.Lts = Lis * parms.gamma;
                    
                    % run simulation
                    if contains(modelname, 'Hill')
                        
                        if strcmp(modelname, 'Hill_alternative')
                            act = parms.actfunc(Cas, parms);
                            
                            Fce = parms.vF_func(parms.vts * parms.gamma, parms) .* act;
                            Fpe = parms.kpe * log(1+exp((parms.Lts - parms.Lce0)));
                            
                            oFi = (Fce + Fpe) * parms.Fscale;
                            
                        else
                            
                            parms.Fse_func = @(dlse, parms) parms.kse0*(exp(parms.kse*dlse)-1);
                            
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
                        end
                        
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
                        %              load('temp4.mat', 'parms')
                        for p = 1:(length(nzi)-1)
                            %                                 disp(p)
                            
                            % simulate
                            if discretized_model
%                                 sol = ode15s(@(t,y,yp) fiber_dynamics_explicit_no_tendon_full(t,y, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, odeopt);
                                sol = ode15s(@(t,y,yp) fiber_dynamics_explicit_discretized(t,y, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, odeopt);
                            else
                                %                                     sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, XP0, []);
                                sol = ode15s(@(t,y) fiber_dynamics_explicit_approximated(t,y, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, []);
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
                                FXB = nan(1, length(sol.x));
                                for iiii = 1:length(sol.x)
                                    n = sol.y(1:end-4,iiii);
                                    xi = parms.xi + L(iiii);
                                    FXB(iiii) = trapz(xi, xi .* n') + trapz(xi, n');
                                end
                                
                                dLce = L - parms.Lce0;
                                Fpe = nan(1, length(sol.x));
                                Fpe((dLce*parms.K) < 10) = parms.kpe * log(1+exp(dLce((dLce*parms.K) < 10)*parms.K))/parms.K;
                                Fpe((dLce*parms.K) >= 10) = parms.kpe * dLce((dLce*parms.K) >= 10);
                                
                                F = FXB + Fpe;
                                
                            else
                                %                                     F = (sol.y(1,:) + sol.y(2,:));
                                F = sol.y(3,:);
                                
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
                        oFi = interp1(tall(ui), Fall(ui), tis) * parms.Fscale; % + parms.Fpe_func(parms.Lts, parms);
                    end
                    %                 toc
                    %%
                    %                         figure(1)
                    %                         plot(tall, Fall)
                    %%
                    % steady state
                    %                         xs = sol.y(:,end);
                    
                    if save_results
                        disp(['Saving ', filename]);
                        cd(output_foldername);
                        save(filename, 'tis','Cas','vis','Lis','oFi','parms','ts')
                        %                             save_model_forces(filename, tis,Cas,vis,Lis,oFi,parms,ts)
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
                        %                             pause
                    end
                end
            end
        end
    end
end


