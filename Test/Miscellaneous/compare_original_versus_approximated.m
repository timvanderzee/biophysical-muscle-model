clear all; close all; 
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

iF = 3;
mcode = [1 2 1];

[output_mainfolder, modelname, ~, ~] = get_folder_and_model(mcode);

parms_version = '_v2';
% parms_version = '';

input_foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iF}];
cd(input_foldername)
load(['parms_',modelname, parms_version, '.mat'])

newparms.kF

%%
% input_foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iF}];
% cd(input_foldername)
% load(['parms_',modelname, '.mat'], 'newparms')
parms = newparms;

pCa = 4.5;
Ca = 10.^(-pCa+6);

AMP = .0383; %.0686;
ISI = .001;

% parms.kpe = 0;

dTt = .0383/.4545; % test stretch (= constant)
dTc = AMP / .4545; % conditioning stretch

% odeopt = [];

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
parms.xi = linspace(-15,15,1000);
% parms.xi = linspace(-25,25,1000);

% parms.f_func = @(xi,f,w)   f/sqrt((2*pi*w^2))*exp(-xi.^2./(2*w^2));
parms.f_func = @(xi,f,w,mu)   f/sqrt((2*pi*w^2))*exp(-(xi-mu).^2./(2*w^2));
parms.g_func = @(xi,k1,k2) k1*exp(k2*xi);

parms.approx = 0;
parms.kF = 1e3;
parms.JF = parms.kF / parms.J1;
parms.ps2 = 0;
parms.b = 5000;
parms.k = 500;

odeopt = [];

for j = 1:2
    tall = [];
    Fall = [];
    Lall = [];

    if j == 1
        x0 = parms.x0';
    elseif j == 2
        x0 = [zeros(length(parms.xi),1); parms.x0(4:end)'];
    end
    
    X0 = x0;
    
    for p = 1:(length(nzi)-1)
        
%           odeopt = odeset('maxstep', 1e-3);
%         if p == 2 && j == 2
%             odeopt = odeset('maxstep', 1e-4);
%         else
%             odeopt = [];
%         end
%         end
        
        disp(p)
        
        xp0 = zeros(size(x0));
        
%         simulate
        if j == 1
            sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, xp0, []);
        else
%              sol = ode15s(@(t,y,yp) fiber_dynamics_explicit_no_tendon(t,y, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, []);
%             parms.k = 500;
%             sol = ode15s(@(t,y,yp) fiber_dynamics_explicit_no_tendon_full(t,y, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, odeopt);
            sol = ode113(@(t,y,yp) fiber_dynamics_explicit_no_tendon_full(t,y, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, odeopt);
        end
        
        t = sol.x;
        X0 = sol.y(:,end);
         
        if j == 2
            n = sol.y(1:length(parms.xi),:);
            L = sol.y(end-3,:);
            
            F = nan(1, length(sol.x));
            for ii = 1:length(sol.x)
                xi = parms.xi + (L(ii) - parms.lce0);
                Qs = trapz(xi(:), [n(:,ii) xi(:).*n(:,ii)]);
                Q0(ii) = Qs(1);
                Q1(ii) = Qs(2);
                F(ii) = sum(Qs);
            end
            

        Lall = [Lall L];
        
        dlse = interp1(parms.ti, parms.Lts, t) - L;
%         F = parms.Fse_func(dlse, parms);
%         [~, id] = unique(tall);
%         

%         F2 = parms.Fse_func(dlse, parms);
        
        else
            F = (sol.y(1,:) + sol.y(2,:));
            
            Q0 = sol.y(1,:);
            Q1 = sol.y(2,:);
            Q2 = sol.y(3,:);
            mu = Q1./Q0; % mean
            q = sqrt(Q2./Q0 - mu.^2);  

            % n
             for ii = 1:length(sol.x)
                 n(:,ii) = Q0(ii) ./ (sqrt(2*pi).*q(ii)) .* exp(-((parms.xi-mu(ii)).^2) ./ (2*q(ii).^2));  % Eq. 52, modified
             end
        end
        
        if p == 2
            
            tids = [2.04 2.05 2.06];
            
            for jj = 1:3
                id = find(t > tids(jj), 1);

                if j == 2
                    xi = parms.xi + (L(id) - parms.lce0);
                else
                    xi = parms.xi;
                end

                figure(100)
                subplot(1,3,jj)
                plot(xi, n(:,id)); hold on
                xlim([-5 5])
            end
            
        end
%         
%         figure(100)
%         subplot(121)
%         plot(t, Q0,'-'); hold on
%         plot(t, Q1, '--')
%         
%         subplot(122)
%         plot(t, F, '-'); hold on
%         plot(t, (Q1+Q0), '--')

        
        tall = [tall t];
        Fall = [Fall F];
        

    end
    
%     if j == 2
%         subplot(311)
%         plot(tall, Lall); hold on
% 
%         subplot(312)
%         plot(tall, Fall); hold on
% 
%         tlin = linspace(0, max(tall), 100);
% 
%         for i = 1:length(parms.xi)
%             ni(i,:) = interp1(tall, n(i,:), tlin);
%         end
% 
%          Li = interp1(tall, L, tlin);
% 
%         figure(10)
%         for i = 1:length(tlin)
%             plot3(tlin(i) * ones(size(parms.xi)), parms.xi, ni(:,i), '-', 'color', [.5 .5 .5]); hold on
%         end
%         ylim([-15 15])
%     end
    
%%
% find unique values
[~, ui] = unique(tall);

% interpolate force
oFi(j,:) = interp1(tall(ui), Fall(ui), tis) * parms.Fscale + parms.Fpe_func(parms.Lts, parms);

end

%%


%%
figure(2)
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


return

%% compare state derivative
parms.xi = linspace(-15,15,2000);

Q0 = 1;
Q1 = .3;
Q2 = .2;

p = Q1./Q0; % mean
q = sqrt(Q2./Q0 - p.^2);  

% n
n = Q0 ./ (sqrt(2*pi)*q) * exp(-((parms.xi-p).^2) / (2*q^2));  % Eq. 52, modified

figure(10)
plot(parms.xi,n)

%%
parms.approx = 0;
L = 0;
Non = .5;
DRX = .5;
R = 1;

% parms.f = 100;
% parms.k11 = 10;
% parms.k21 = 10;
% 
% parms.ps2 = 1;
% parms.dLcrit = 10;
% 
% parms.kse = 1;
% parms.kse0 = 1;

y1 = [Q0 Q1 Q2 L Non DRX R]';
[dx1] = fiber_dynamics_explicit_no_tendon(0,y1, parms);

y2 = [n L Non DRX R]';
[dx2] = fiber_dynamics_explicit_no_tendon_full(0,y2, parms);

%%
clc
Ld = dx1(end-3);
Qdot0 = trapz(parms.xi, dx2(1:length(parms.xi)));
Qdot1 = trapz(parms.xi, parms.xi.*dx2(1:length(parms.xi))') + 1 * Ld .* Q0;
Qdot2 = trapz(parms.xi, parms.xi.^2.*dx2(1:length(parms.xi))') + 2 * Ld .* Q1;

disp([Qdot0 Qdot1 Qdot2])
disp(dx1(1:3))


%%




