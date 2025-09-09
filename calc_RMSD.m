clear all; close all; clc
usernames = {'timvd','u0167448'};

for i = 1:length(usernames)
    docfolder =  ['C:\Users\', usernames{i}, '\Documents'];
    
    if isfolder(docfolder)
        mainfolder = docfolder;
    end
end

if contains(mainfolder, 'timvd')
    githubfolder = mainfolder;
else
    githubfolder = [mainfolder, '\GitHub'];
end

addpath(genpath([githubfolder, '\muscle-thixotropy']))
addpath(genpath([mainfolder, '\casadi-3.7.1-windows64-matlab2018b']))
addpath(genpath([githubfolder, '\biophysical-muscle-model']))

mcodes = [1 1 1; 1 1 1; 1 1 3; 2 1 1];
      
for kk = 1:size(mcodes,1)
    mcode = mcodes(kk,:);
    
    [output_mainfolder, filenames{kk}, opt_types{kk}, ~] = get_folder_and_model(mcode);
end

%% load data and parameters
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

cd([output_mainfolder{2},'\data'])

k = 7;
load([fibers{k},'_cor_new.mat'],'data');

% chosen ISIs, AMPs and pCas
ISIs = repmat([.001 .100 .316 1], 4,1);

AMPs = repmat([0 .0038 .0121 .0383]', 1, 4); 

pCas = 6.1 * ones(size(AMPs));

vs = {'\', '\'};

for kk = 1:size(mcodes,1)
    disp(filenames{kk})
    output_folder = [opt_types{kk},'\normalized\with_PE_optimized\2_trials'];

    output_dir = [output_mainfolder{1}, '\', filenames{kk},vs{1}, output_folder];
    cd(output_dir)

    load([filenames{kk},'_F', num2str(k),'_best.mat'],'parms','exitflag','fopt','C0','Cbounds','model','P0','P')
    C = p_to_c(P, Cbounds);
    parms = C_to_parms(C, parms, parms.optvars);
    parms = calc_dependent_parms(parms);  

    parms.act = 1;
    parms.Noverlap = 1;
    parms.vF_func =  @(vcerel,parms)parms.e(1)*log((parms.e(2)*vcerel./parms.vmax+parms.e(3))+sqrt((parms.e(2)*vcerel./parms.vmax+parms.e(3)).^2+1))+parms.e(4);

    Parms{kk} = parms;
    
    if kk == 1
        cd([mainfolder, '\GitHub\biophysical-muscle-model\Parameters'])
        load('parms_v5.mat')
        Parms{kk} = sparms(k);
    end
end

%% evaluate
odeopt = odeset('maxstep', 1e-2);
gamma = 108.3333; % length scaling

RMSD = nan(size(ISIs,1), size(ISIs,2), size(mcodes,1));

for j = 1:size(ISIs,1)

    [texp, Lexp, Fexp, Tsrel] = get_data(data, ISIs(j,:), AMPs(j,:), pCas(j,:));

    for i = 1:size(ISIs,2)
        
        % evaluate fit
        Ca = 10.^(-pCas(j,i)+6);

        dTt = .0383/.4545; % test stretch (= constant)
        dTc = AMPs(j,i) / .4545; % conditioning stretch
        ISI = ISIs(j,i);
        tiso = 3;

        [tis, Cas, Lis, vis, ts] = create_input(tiso, dTt, dTc, ISI, Ca, 2000);
        
        % center around 2nd stretch
        t = tis + 3*dTt - tiso;
        
        for kk = 1:size(mcodes,1)
        
            parms = Parms{kk};

            x0 = parms.x0(2:end)';
            xp0 = zeros(size(x0));

            parms.ti = tis;
            parms.vts = vis;
            parms.Cas = Cas;
            parms.Lts = Lis;

            % run simulation
            tic
            if contains(filenames{kk}, 'Hill')
                sol = ode15i(@(t,y,yp) hill_type_implicit(t,y,yp, parms), [0 max(tis)], x0, xp0, odeopt);
            else
                sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 max(tis)], x0, xp0, odeopt);
            end
            toc

            % update x0
            x0 = sol.y(:,end);

            Liss = Lis * gamma;
            ot = sol.x;

            if contains(filenames{kk}, 'Hill')
                lce = sol.y(1,:);

                % elastic elements
                Lse = Liss - interp1(sol.x, lce, tis);    
                oFi = parms.Fse_func(Lse, parms) * parms.Fscale + parms.Fpe_func(Liss, parms);

            else
                F = sol.y(1,:) + sol.y(2,:);

                oF = F * parms.Fscale;

                oFi = interp1(ot, oF, tis) + parms.Fpe_func(Liss, parms);
            end

            % compute RMSD
            id = t < .15 & t > (-ISI - 2 * dTc - .1);
            oFii = interp1(t(id), oFi(id), texp(:,i));
            
            RMSDs = sqrt((oFii - Fexp(:,i)).^2) * 100;
            
            RMSD(j,i,kk) = sqrt(mean((oFii - Fexp(:,i)).^2, 'omitnan'));
%             axis([-.35 .25 0 20])
        end
    end
end

%% 
close all
figure(1)

for kk = 1:size(mcodes,1)
nexttile
surf(ISIs, AMPs, RMSD(:,:,kk))
xlabel('ISI')
ylabel('AMP')

set(gca, 'XScale', 'log')
% zlim([0 2])
end