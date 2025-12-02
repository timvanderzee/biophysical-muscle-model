clear all; close all ;clc
[username, githubfolder] = get_paths();
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

iFs = 7;
parms_version = '_v3';

% load parameters
mcode = [1 1 1];
[output_mainfolder, modelname, ~, ~] = get_folder_and_model(mcode);

input_foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iFs}];
cd(input_foldername)
load(['parms_',modelname, parms_version, '.mat'], 'newparms')
parms = update_parms(newparms);


pCa     = 6.1;
Ca      = 10.^(-pCa+6);
AMP     = .0383;
ISI     = .001;
dTt     = .0383/.4545; % test stretch (= constant)
dTc     = AMP / .4545; % conditioning stretch
tiso    = dTt*3+dTc*2+ISI + 2;
dt      = .001; % gives 10 points in SRS zone
N       = round(tiso / dt);
[tis, Cas, Lis, vis, ts, Ts] = create_input(tiso, dTt, dTc, ISI, Ca, N);

%% simulate entire protocol
% discretization parameters
bw = .03; % standard bin width
Ns = [logspace(1,3,7) 1e4]; % bin#

% experiment
tlin = linspace(0,Ts(end-1),1000);
dt = 1e-3;

parms.Cas = mean(Cas);

% pre-allocate
tsim = nan(length(Ns),2,2);
Fi = nan(length(tlin), length(Ns));

% tvec = 0:dt:max(tlin);
% N = length(tvec);
% vts = interp1(tis, vis, tvec);

for f = 1:3 % methods
    if f == 1
        parms.method_of_characteristics = 0;
    else
        parms.method_of_characteristics = 1;
    end

    for j = 1:2
        for id = 1:length(Ns)
            
            if f < 3
                if j == 1
                    parms.xi = linspace(-15,15, round(Ns(id)));
                else
                    parms.xi = linspace(-bw*round(Ns(id)/2),bw*round(Ns(id)/2), round(Ns(id)));
                end

                n0 = zeros(size(parms.xi));
                X0 = [n0'; parms.x0(4:end)'];

            else
                X0 = parms.x0'; 
            end
            
            t2 = 0;
            F = 0;
            y = X0;
            
            disp(id)
           
            tic
            ii = 1;
            
            while t2(ii) <Ts(end)
                ii = ii+1;
                
                parms.F = F(ii-1);
%                 parms.vts = vts(ii);
                parms.vts = interp1(tis, vis, t2(ii-1));
                
                if f < 3
                    [yp, F(ii), Q0] = fiber_dynamics_explicit_no_tendon_full(t2(ii-1),y(:,ii-1), parms);
                else
                    [yp, F(ii), Q0] = fiber_dynamics_explicit_no_tendon(t2(ii-1),y(:,ii-1), parms);
                end
                
%                 if t2(ii-1) < .01
%                     dt = 3e-3;
%                 else
%                     dt = 3e-3;
%                 end
                
                y(:,ii) = y(:,ii-1) + yp * dt;
                t2(ii) = t2(ii-1) + dt;
                  
                if f == 1 % shift XB distribution
                    % for n, also need to shift
                    dx = yp(end-3) * dt;
                    xi = parms.xi - dx;
                    n = y(1:length(parms.xi),ii);
                    
                    nshift = interp1(parms.xi(:), n, xi(:),'linear',0);
                    nshift(isnan(nshift)) = 0;
                    nshift(nshift<0) = 0;
                    
                    y(1:length(parms.xi),ii) = nshift;
   
%                     Qs = trapz(xi(:), [nshift(:) xi(:).*nshift(:)]);
%                     F(ii) = sum(Qs);
                    
                     F(ii) = F(ii) + dx * Q0;
                end
            end
        
            tsim(id,j,f) = toc;

            % compute force
            Fi(:,id,j,f) = interp1(t2, F, tlin);

        end
    end
end

N = ii;

%% forces
if ishandle(2), close(2); end; figure(2)
% 
% % subplot(221)
plot(tlin,Fi(:,end,1,2)); hold on
plot(tlin,Fi(:,end,2,2),'--')
plot(tlin,Fi(:,end,1,3),':')

% plot(tlin, Fapi,':')

% subplot(222)
% semilogy(tlin, abs(Fi - Fi(:,end)));

%%
close all
titles = {'Fixed strain range', 'Fixed bin width'};

figure(3)

color = get(gca,'colororder');

Ntot = 1 * ones(size(tsim));

for j = 1:2
    for f = 1:3
        
        subplot(3,2,j)
        semilogx(Ns, tsim(:,j,f) ./ Ntot(:,j,f),'ko--', 'markerfacecolor', color(f,:)); hold on
        box off
        xlabel('# bins')
        ylabel('Cost (s)')
        
        title(titles{j})
        subtitle('Cost')
        
        xlim([Ns(1) Ns(end-1)])
        
        ylim([0 max(tsim(1:end-1,:,:),[],'all')/mean(Ntot(:), 'omitnan')])
        %
%         if f == 2
%             yline(tapi/Napi, 'k--', 'linewidth', 1)
%         end
        
        subplot(3,2,j+2)
        semilogx(Ns, mean(abs(Fi(:,:,j,f) - Fi(:,end,j,f)), 'omitnan'),'ko--', 'markerfacecolor', color(f,:)); hold on
        box off
        subtitle('Error')
        xlabel('# bins')
        ylabel('Error (-)')
        
        xlim([Ns(1) Ns(end-1)])
        ylim([1e-10 .4])
        
%         if f == 2
%             yline(mean(abs(Fapi(:) - Fi(:,end,j,f)), 'omitnan'),'k--', 'linewidth', 1)
%         end
        
        subplot(3,2,j+4)
        plot(tsim(1:end-1,j,f) ./ Ntot(1:end-1,j,f), mean(abs(Fi(:,1:end-1,j,f) - Fi(:,end,j,f)), 'omitnan'),'ko--', 'markerfacecolor', color(f,:)); hold on
        box off
        subtitle('Error versus cost')
        xlabel('Cost (s)')
        ylabel('Error (-)')
        
%         if f == 2
%             plot(tapi/Napi, mean(abs(Fapi(:) - Fi(:,end,j,f)), 'omitnan'), 'ko', 'markerfacecolor', [0 0 0])
%             yline(mean(abs(Fapi(:) - Fi(:,end,j,f)), 'omitnan'),'k--', 'linewidth', 1)
%             xline(tapi/Napi, 'k--', 'linewidth', 1)
%         end
        
        ylim([0 .01])
        xlim([0 max(tsim(1:end-1,:,:),[],'all')/mean(Ntot(:), 'omitnan')])
        %         xlim([Ns(1) Ns(end-1)])
        
    end
end

%%

subplot(321)
legend('Tradit.', 'Charac.', 'Approx.', 'location', 'best')
legend boxon


% subplot(326)
% legend('Tradit.', 'Charac.', 'Approx.', 'location', 'best')
% legend boxon


%%
% interval needs to have finite duration
aTs = [0; Ts];
nzi = find(diff(aTs) > 0);
vs = [0 .4545 -.4545 0 .4545 0 0];
X0 = parms.x0'; 
Y = [];
t1 = [];

% for p = 1:(length(nzi)-1)
% 
%     parms.vts = vs(nzi(p));
% 
%     sol = ode45(@(t,y,yp) fiber_dynamics_explicit_no_tendon(t,y, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, []);
% 
%     Y = [Y sol.y(:,1:end-1)];
%     t1 = [t1 sol.x(1:end-1)];
% 
%     X0 = sol.y(:,end);
% 
% %     Nit(p) = sol.stats.nfevals;
% 
% end
% Fap = Y(1,:) + Y(2,:);

parms.vts = vis;
parms.ti = tis;
sol = ode45(@(t,y,yp) fiber_dynamics_explicit_no_tendon(t,y, parms), [0 aTs(end)], X0, []);

t1 = sol.x;
Fap = sol.y(1,:) + sol.y(2,:);


%%
figure(1)

subplot(211)
semilogy(t1(1:end-1), diff(t1))
title('Choosen time step by ode45')
xlabel('Time (s)')
ylabel('Time step (s)')
box off
yline(6e-3,'k--')

subplot(212);
plot(t1, Fap)
title('Force')
xlabel('Time (s)')
box off
ylabel('Force')