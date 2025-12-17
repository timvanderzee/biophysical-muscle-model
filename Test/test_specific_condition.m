clear all; close all; clc

[username, githubfolder] = get_paths();

mcodes = [1 2 1; 1 2 1];
discretized_model = [0 1];
FLs = [0 1];

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

iFs = 6;
pCa = 9;
Ca = 10.^(-pCa+6);
AMP = .0383;
ISI = .001;
parms_version = '_v4';

figure(1)
for iF = iFs
    nexttile
    for i = 1:length(discretized_model)

        % load parameters
        mcode = mcodes(i,:);
        [output_mainfolder, modelname, ~, ~] = get_folder_and_model(mcode);
        
        input_foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iF}];
        cd(input_foldername)
        load(['parms_',modelname, parms_version, '.mat'], 'newparms')
        parms = update_parms(newparms);
%         parms.FL_overlap = FLs(i);
  
      % find L0
        costfunc = @(L, Fce, parms) (Fce + (parms.kpe*(-L-parms.Lce0)) - (parms.kse0*(exp(parms.kse*L)-1))) .^2;
        dlse = fminsearch(@(L) costfunc(L, 0, parms), 0);

        lce0 = -dlse;
        Fse0 = parms.Fse_func(dlse, parms);
        
%         load('temp4.mat')
        
        if contains(modelname, 'Hill')
            x0 = 0;
        elseif discretized_model(i)
            
            n0 = zeros(size(parms.xi));
            x0 = [n0'; parms.x0(4:end)'];

  
            x0(end-3) = lce0;
%             x0(end-3) = 0;
            
        else
            x0 = parms.x0';
%             x0(3) = Fse0;
            
            Q0 = .1e-6;
            p0 = -1;
            q0 = .1;
            [Q00, Q20, lce0, Q10, Fse0] = find_steady_state(Q0, p0, q0, parms, 'regular');
            x0 = [Q00 Q20 Fse0 0 .1 0];
            
        end
        
%         load('temp4.mat')
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
                    
                    %                 sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, xp0, []);
                    sol = ode15s(@(t,y) fiber_dynamics_explicit_length_v2(t,y, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, odeopt);
                    
                    % get force
                    t = sol.x;
                    F = sol.y(3,:);
                    L = interp1(parms.ti, parms.Lts, t) - parms.Lse_func(F, parms);
                    
                end
                
%                 Fmax(p) = max(F);
                X0 = sol.y(:,end);

                tall = [tall t];
                Fall = [Fall F];
                Lall = [Lall L];
                
            end
            
            % find unique values
            [~, ui] = unique(tall);
            
            % interpolate force
            oFi = interp1(tall(ui), Fall(ui), tis) * parms.Fscale; % + parms.Fpe_func(parms.Lts, parms);
        end
        
        
        figure(1)
        subplot(211)
        plot(tall, Fall * parms.Fscale, 'linewidth', 2); hold on
        box off
        xlabel('Time (s)')
        ylabel('Force (-)')
        yline(1,'k--')
        
        subplot(212)
        plot(tall, Lall, 'linewidth', 2); hold on
        box off
        xlabel('Time (s)')
        ylabel('Length (-)')
%         xlim([4.5 max(tall)])
    end
end


