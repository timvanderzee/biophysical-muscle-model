clear all; close all; clc

[username, githubfolder] = get_paths();

mcodes = [1 2 1; 1 2 1; 1 2 1];
discretized_model = [0 0 1];

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

iFs = 2;
pCa = 4.5;
Ca = 10.^(-pCa+6);
AMP = .0383;
ISI = .001;
parms_version = '_v2';

figure(1)
for iF = iFs
    nexttile
    for i = 1:size(mcodes,1)

        % load parameters
        mcode = mcodes(i,:);
        [output_mainfolder, modelname, ~, ~] = get_folder_and_model(mcode);
        
        input_foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iF}];
        cd(input_foldername)
        load(['parms_',modelname, parms_version, '.mat'], 'newparms')
        
        parms = newparms;
        
        parms.xi = linspace(-15,15,1000);
        parms.f_func = @(xi,f,w,mu)   f/sqrt((2*pi*w^2))*exp(-(xi-mu).^2./(2*w^2));
        parms.g_func = @(xi,k1,k2) k1*exp(k2*xi);
        parms.approx = 1;
        
%         if i == 1 || i == 3
            parms.k = 1e3;
            parms.b = 5e3;
            parms.dLcrit = 0.7;
            parms.ps2 = parms.dLcrit-parms.w;
%         end

        if i == 2
            parms.approx = 1;
        else
            parms.approx = 0;
        end
            
        if contains(modelname, 'Hill')
            x0 = 0;
        elseif discretized_model(i)
            
            n0 = zeros(size(parms.xi));
            x0 = [n0'; parms.x0(4:end)'];
            
        else
            x0 = parms.x0';
        end
        
        xp0 = zeros(size(x0));
        X0 = x0;
       
        dTt = .0383/.4545; % test stretch (= constant)
        dTc = AMP / .4545; % conditioning stretch
        
        odeopt = odeset('maxstep', 1e-3);
        odeopt = [];
        
        tiso = dTt*3+dTc*2+ISI + 5;
        dt = .001; % gives 10 points in SRS zone
        N = round(tiso / dt);

        [tis, Cas, Lis, vis, ts, Ts] = create_input(tiso, dTt, dTc, ISI, Ca, N);
        
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
            oFi = parms.Fse_func(Lse, parms) * parms.Fscale + parms.Fpe_func(parms.Lts, parms);
            
        else
            
            aTs = [0; Ts];
            X0 = x0;
            vts = [0 .4545 -.4545 0 .4545 0 0];
            
            % splitting it up makes things much faster
            tall = [];
            Fall = [];
            Lall = [];
            
            % interval needs to have finite duration
            nzi = find(diff(aTs) > 0);
            
            for p = 1:(length(nzi)-1)
                disp(p)
                
                % simulate
                if discretized_model(i)
                    sol = ode15s(@(t,y,yp) fiber_dynamics_explicit_no_tendon_full(t,y, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, odeopt);
                    t = sol.x;
                    
                    L = sol.y(end-3,:);
                    
                    F = nan(1, length(sol.x));
                    for iiii = 1:length(sol.x)
                        n = sol.y(1:end-4,iiii);
                        xi = parms.xi + L(iiii);
                        F(iiii) = trapz(xi, xi .* n') + trapz(xi, n');
                    end
                    
                else
                    
                    %                 sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, xp0, []);
                    sol = ode15s(@(t,y) fiber_dynamics_explicit_no_tendon(t,y, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, odeopt);
                    
                    % get force
                    t = sol.x;
                    F = (sol.y(1,:) + sol.y(2,:));
                    L = interp1(parms.ti, parms.Lts, t) - parms.Lse_func(F, parms);
                    
                end
                
                  X0 = sol.y(:,end);

                tall = [tall t];
                Fall = [Fall F];
                Lall = [Lall L];
                
            end
            
            % find unique values
            [~, ui] = unique(tall);
            
            % interpolate force
            oFi = interp1(tall(ui), Fall(ui), tis) * parms.Fscale + parms.Fpe_func(parms.Lts, parms);
        end
        
        
        figure(1)
        plot(tis, oFi, 'linewidth', 2); hold on
        box off
        xlabel('Time (s)')
        ylabel('Force (-)')
        
%         xlim([4.5 max(tall)])
    end
end


