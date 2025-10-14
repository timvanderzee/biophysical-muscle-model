clear all; close all; clc
[username, githubfolder] = get_paths();

mcodes = [1 2 1; 2 1 1];
iFs = 6;
pCas = [4.5 6.1 6.2 6.3 6.4 6.6 9];
pCas = 6.2;

Ca = 10.^(-pCas+6);

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

visualize = 0;

%% create input
tis = linspace(0,3,1000);
f = 2;

AMP = .04;
Lis = AMP - AMP * cos(2*pi*f*tis);
vis = 2*pi*f*AMP * sin(2*pi*f*tis);

Lis(tis<1) = 0;
vis(tis<1) = 0;


%% simulate


for iii = 1:size(mcodes,1)
    
    % load parameters
    mcode = mcodes(iii,:);
    [output_mainfolder, filename, ~, ~] = get_folder_and_model(mcode);
    
    % disp(filename)
    
    cd([githubfolder, '\biophysical-muscle-model\Parameters'])
    load(['parms_',filename,'.mat'], 'pparms')
    
    gamma = 108.3333; % length scaling
    
    %% step 1: force - pCa
    parms = pparms(1);
    
    if contains(filename, 'Hill')
        x0 = 0;
    else
        x0 = parms.x0';
        
    end
    
    xp0 = zeros(size(x0));
    
    
    %% evaluate
    for iF = iFs
        
        parms = pparms(iF);
        parms.K = 100;
        
        
        for i = 1:length(Ca)
            
            X0 = x0;
            
            
            parms.ti = tis;
            parms.vts = vis;
            parms.Cas = Ca(i);
            parms.Lts = Lis;
            Liss = Lis * gamma;
            
            % run simulation
            if contains(filename, 'Hill')
                % simulate
                sol = ode15i(@(t,y,yp) hill_type_implicit(t,y,yp, parms), [0 max(tis)], X0, xp0, []);
                
                % get SE length
                Lse = Liss - interp1(sol.x, sol.y(1,:), tis);
                
                % get force
                oFi = parms.Fse_func(Lse, parms) * parms.Fscale + parms.Fpe_func(Liss, parms);
                
            else
                
                % simulate
                sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 max(tis)], X0, xp0, []);

                X0 = sol.y(:,end);

                % get force
                t = sol.x;
                F = (sol.y(1,:) + sol.y(2,:));
                    
                % interpolate force
                oFi = interp1(t, F, tis) * parms.Fscale + parms.Fpe_func(Liss, parms);
            end
            
            figure(1)
            subplot(211)
            plot(tis, Lis,'linewidth',2, 'color', [.5 .5 .5]); hold on
            
            subplot(212)            
            plot(tis, oFi,'linewidth',2); hold on
        end
    end
end

%% make nice
ylabels = {'Length (L_0)', 'Force (F_0)'};

ymax = [.08 .6];

figure(1)

for i = 1:2
    subplot(2,1,i)
    box off
    xlim([.8 2.5])
    plot([1 1], [0 ymax(i)],'k--')
    plot([1.5 1.5], [0 ymax(i)],'k--')
    plot([2 2], [0 ymax(i)],'k--')
    plot([2.5 2.5], [0 ymax(i)],'k--')

    xlabel('Time (s)')
    ylabel(ylabels{i})
end

subplot(212)
legend('Biophysical','Hill-type','location','best')
legend boxoff

    

    
    
    
    
    
    
    
    
    
    