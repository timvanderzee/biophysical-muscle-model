clear all; close all; clc
save_results = 1;

[username, githubfolder] = get_paths();

mcodes = [1 1 1; 1 2 1];

iFs = [2,3,5,6,7,8,11];
AMPs = [0    0.0012    0.0038    0.0121    0.0216    0.0288    0.0383    0.0532    0.0682];
ISIs = [ 0.0010    0.0100    0.0500    0.1000    0.2000    0.3160    0.5000    1.0000    3.1600   10.0000];
pCas = [4.5 6.1 6.2 6.3 6.4 6.6 9];
Ca = 10.^(-pCas+6);
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

visualize = 0;

iF = 7;
mcode= [1 1 1];

[output_mainfolder, modelname, ~, ~] = get_folder_and_model(mcode);

input_foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iF}];
cd(input_foldername)
load(['parms_',modelname, '.mat'], 'newparms')
parms = newparms;

pCa = 9;
Ca = 10.^(-pCa+6);

AMP = .0383;
ISI = 1;

parms.kpe = 0;

dTt = .0383/.4545; % test stretch (= constant)
dTc = AMP / .4545; % conditioning stretch

odeopt = [];

tiso = dTt*3+dTc*2+ISI + 2;
dt = .001; % gives 10 points in SRS zone
N = round(tiso / dt);

[tis, Cas, Lis, vis, ts, Ts] = create_input(tiso, dTt, dTc, ISI, Ca, N);

parms.ti = tis;
parms.vts = vis;
parms.Cas = mean(Cas);
parms.Lts = Lis * parms.gamma;

% run simulation
aTs = [0; Ts];
vts = [0 .4545 -.4545 0 .4545 0 0];

% splitting it up makes things much faster


%% test models
% interval needs to have finite duration
nzi = find(diff(aTs) > 0);
parms.xi = linspace(-15,15,500);
parms.f_func = @(xi,f,w)   f/sqrt((2*pi*w^2))*exp(-xi.^2./(2*w^2));
parms.g_func = @(xi,k1,k2) k1*exp(k2*xi);

for j = 1:2
    tall = [];
    Fall = [];
    
    if j == 1
        x0 = parms.x0';
    elseif j == 2
        x0 = [zeros(length(parms.xi),1); parms.x0(4:end-1)'];
    end
    
    X0 = x0;
    
    for p = 1:(length(nzi)-1)
        
        disp(p)
        
        xp0 = zeros(size(x0));
        
        % simulate
        if j == 1
            sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, xp0, []);
        else
            sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon_full(t,y,yp, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, xp0, []);
        end
        
        t = sol.x;
        X0 = sol.y(:,end);
         
        if j == 2
            n = sol.y(1:length(parms.xi),:);
            L = sol.y(end-2,:);
            
            F = nan(1, length(sol.x));
            for ii = 1:length(sol.x)
                xi = parms.xi + (L(ii) - parms.lce0);
                Qs = trapz(xi(:), [n(:,ii) xi(:).*n(:,ii)]);
                F(ii) = sum(Qs);
            end
            
        else
            F = (sol.y(1,:) + sol.y(2,:));
        end
        
        tall = [tall t];
        Fall = [Fall F];
        
        % simulate
        %     solf = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon_full(t,y,yp, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, xp0, []);
        
        
    end
    
    plot(tall, Fall); hold on


%%
% find unique values
[~, ui] = unique(tall);

% interpolate force
oFi(j,:) = interp1(tall(ui), Fall(ui), tis) * parms.Fscale + parms.Fpe_func(parms.Lts, parms);

end

%%

figure(1)
subplot(2,1,1)
plot(tis, parms.Lts, 'linewidth', 2); hold on
box off
xlabel('Time (s)')
ylabel('Length (-)')

subplot(2,1,2)
plot(tis, oFi, 'linewidth', 2); hold on
box off
xlabel('Time (s)')
ylabel('Force (-)')

