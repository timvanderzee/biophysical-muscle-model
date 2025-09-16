clear all; close all; clc
[username, githubfolder] = get_paths();

% load parameters
mcode = [2 1 1];
[output_mainfolder, filename, ~, ~] = get_folder_and_model(mcode);

cd([githubfolder, '\biophysical-muscle-model\Parameters'])
load(['parms_',filename,'.mat'], 'pparms')

gamma = 108.3333; % length scaling

%% step 1: force - pCa
pCas = [4.5 6.1 6.2 6.3 6.4 6.6 9];
Ca = 10.^(-pCas+6);

parms = pparms(1);

if contains(filename, 'Hill')
    x0 = 0;
else
    x0 = [parms.x0(2:end)'; 0];
end

xp0 = zeros(size(x0));


%% evaluate
iFs = [1 2 3, 5,6, 7, 8, 10, 11];

AMPs = [0 12 38 121 216 288 383 682]/10000;
% ISIs = [1 10 100 316 1000]/1000;
ISIs = [1 10 100 316 1000 3160 10000]/1000;

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

visualize = 1;

for iF = 6
    
    parms = pparms(iF);
    
    for i = 1:length(Ca)
          
        % initial condition
%         x0 = xs(:,i,iF);       
%         xp0 = zeros(size(x0));
        
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
                
                tiso = dTt*3+dTc*2+ISI + 2;
                
                dt = .001; % gives 10 points in SRS zone
                N = round(tiso / dt);
                
                [tis, Cas, Lis, vis, ts] = create_input(tiso, dTt, dTc, ISI, Ca(i), N);
                
                parms.ti = tis;
                parms.vts = vis;
                parms.Cas = mean(Cas);
                parms.Lts = Lis;
                Liss = Lis * gamma;
                
                % run simulation
                if contains(filename, 'Hill')
                    % simulate
                    sol = ode15i(@(t,y,yp) hill_type_implicit(t,y,yp, parms), [0 max(tis)], x0, xp0, odeopt);
                    
                    % get SE length
                    Lse = Liss - interp1(sol.x, sol.y(1,:), tis);
                    
                    % get force
                    oFi = parms.Fse_func(Lse, parms) * parms.Fscale + parms.Fpe_func(Liss, parms);
                    
                else
                    % simulate
                    sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 max(tis)], x0, xp0, odeopt);
                    
                    % get force
                    F = (sol.y(1,:) + sol.y(2,:));
                    
                    % interpolate force
                    oFi = interp1(sol.x, F, tis) * parms.Fscale + parms.Fpe_func(Liss, parms);
                end
                
                cd([output_mainfolder{2}])
                
                if ~isfolder([filename,'\',fibers{iF}, '\pCa=',num2str(pCas(i)*10)])
                    mkdir([filename,'\',fibers{iF}, '\pCa=',num2str(pCas(i)*10)])
                end
                
                cd([filename,'\',fibers{iF}, '\pCa=',num2str(pCas(i)*10)])
                disp([filename,'\',fibers{iF}, '\pCa=',num2str(pCas(i)*10),'\', fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'])
                save([fibers{iF},'_AMP=',num2str(AMP*10000),'_ISI=',num2str(ISI*1000),'.mat'], ...
                    'tis','Cas','vis','Lis','oFi','parms','ts')
                
                if visualize
                    figure(1)
                    subplot(411)
                    plot(tis, Cas, 'b', 'linewidth',1);
                    
                    subplot(412)
                    plot(tis, vis,'b',  'linewidth',1);
                    
                    subplot(413)
                    plot(tis, Lis,'b',  'linewidth',1);
                    
                    subplot(414);
                    plot(tis, oFi,'b')
                    
                    drawnow
                end
            end
        end
    end
end

