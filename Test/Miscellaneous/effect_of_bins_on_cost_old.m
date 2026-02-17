clear all; close all ;clc
[username, githubfolder] = get_paths();

discretized_model = 1;
FLs = [0 1];

fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};

iFs = 6;
pCa = 6.1;
Ca = 10.^(-pCa+6);
AMP = 0;
ISI = .001;
parms_version = '_v3';

% load parameters
mcode = [1 1 1];
[output_mainfolder, modelname, ~, ~] = get_folder_and_model(mcode);

input_foldername = [githubfolder, '\biophysical-muscle-model\Parameters\',fibers{iFs}];
cd(input_foldername)
load(['parms_',modelname, parms_version, '.mat'], 'newparms')
parms = update_parms(newparms);
%     parms.FL_overlap = FLs(i);

n0 = zeros(size(parms.xi));
x0 = [n0'; parms.x0(4:end)'];
xp0 = zeros(size(x0));
X0 = x0;

dTt = .0383/.4545; % test stretch (= constant)
dTc = AMP / .4545; % conditioning stretch

odeopt = odeset('maxstep', 1e-3);

tiso = dTt*3+dTc*2+ISI + 2;
dt = .001; % gives 10 points in SRS zone
N = round(tiso / dt);

[tis, Cas, Lis, vis, ts, Ts] = create_input(tiso, dTt, dTc, ISI, Ca, N);

%% simulate
bw = mean(diff(linspace(-15,15,500)));

clc

parms.vts = 0;
parms.Cas = 1;
parms.method_of_characteristics = 1;

K = 1e4;

Ns = [51 251 501 1001 5001];
N = length(Ns);

tfval = nan(N,K);
tint = nan(N,K);

close all
figure(1)

for j = 1:2
    for k = 1:K
        for id = 1:N
            
            %         disp(id)
            
            if j == 1
                parms.xi = linspace(-3,3, round(Ns(id)));
            else
                parms.xi = linspace(-bw*round(Ns(id)/2),bw*round(Ns(id)/2), round(Ns(id)));
            end
            
            n0 = zeros(size(parms.xi));
            X0 = [n0'; zeros(4,1)];
            
            tic
            yp = fiber_dynamics_explicit_no_tendon_full(0,X0, parms);
            tfval(id,k,j) = toc;
            
            Nx(id) = length(X0);
            
            % interpolating
            tic
            nshift = interp1(parms.xi(:), n0, parms.xi(:)+.1);
            %         Qs = trapz(parms.xi(:), n0(:).*parms.xi(:));
            tint(id,k,j) = toc;
            
        end
    end
end


%% plot
tfval(:,1,:) = nan;
tint(:,1,:) = nan;

figure(1)
titles = {'Fixed range', 'Fixed bin width'};

for j = 1:2
    
    ttot = tfval + tint;
    
    % fit evaluations only
    
    
    % semilogx([50 50], [0 1e-3],'r:', 'linewidth',1); hold on
    % semilogx([500 500], [0 1e-3],'g:', 'linewidth',1)
    color = get(gca, 'colororder');
    
    % plot(Ns, tfval, 'color', [.5 .5 .5 .9]); hold on
    Nall = logspace(1, 4, 1000);
    
    subplot(1,2,j)
    errorbar(Ns, mean(tfval(:,:,j),2, 'omitnan'), std(tfval(:,:,j),1,2,  'omitnan'),'o', 'color',[.5 .5 .5], 'markerfacecolor', [.5 .5 .5]); hold on
    errorbar(Ns, mean(ttot(:,:,j),2, 'omitnan'), std(ttot(:,:,j),1,2,  'omitnan'),'o', 'color', color(1,:), 'markerfacecolor', color(1,:))
    
    set(gca,'XScale','log')
    
    p = polyfit(Ns, mean(tfval(:,:,j),2, 'omitnan'), 1);
    semilogx(Nall, polyval(p, Nall), 'k--', 'linewidth',1); hold on
    
    p = polyfit(Ns, mean(ttot(:,:,j),2, 'omitnan'), 1);
    semilogx(Nall, polyval(p, Nall), 'k--', 'linewidth',1); hold on
    
    errorbar(Ns, mean(tfval(:,:,j),2, 'omitnan'), std(tfval(:,:,j),1,2,  'omitnan'),'o', 'color',[.5 .5 .5], 'markerfacecolor', [.5 .5 .5]); hold on
    errorbar(Ns, mean(ttot(:,:,j),2, 'omitnan'), std(ttot(:,:,j),1,2,  'omitnan'),'o', 'color', color(1,:), 'markerfacecolor', color(1,:))
    
    
    legend('Excl. interpolation', 'Incl. interpolation', 'location', 'best')
    legend boxoff
    grid on
    
    box off
    xlabel('# bins')
    ylabel('Cost per function evaluation (s)')
    xlim([20 1e4])
    ylim([0 1e-3])
    
    title(titles{j})
end

%% simulate experiment
bw = mean(diff(linspace(-15,15,500)));

% interval needs to have finite duration
aTs = [0; Ts];
nzi = find(diff(aTs) > 0);
vs = [0 .4545 -.4545 0 .4545 0 0];

tlin = linspace(0,Ts(end-1),1000);

parms.Cas = mean(Cas);

Ns = logspace(1,3,8);
tsim = nan(length(Ns),2,2);
Ntot = nan(length(Ns),2,2);
Fi = nan(length(tlin), length(Ns));
Nit = nan(length(nzi)-1, length(Ns), 2);

dt = 1e-3;
tvec = 0:dt:max(tlin);
N = length(tvec);
vts = interp1(tis, vis, tvec);

for f = 1 % methods
    for j = 1
        for id = length(Ns)
            if j == 1
                parms.xi = linspace(-15,15, round(Ns(id)));
            else
                parms.xi = linspace(-bw*round(Ns(id)/2),bw*round(Ns(id)/2), round(Ns(id)));
            end
            
            n0 = zeros(size(parms.xi));
            X0 = [n0'; parms.x0(4:end)'];
            
            tic
            
            t2 = 0;
            F2 = 0;
            y = X0;
            
            disp(id)
            
            if f == 1
                % method of characteristics
                parms.method_of_characteristics = 0;
                for ii = 2:N
                    
                    parms.F = F2(ii-1);
                    parms.vts = vts(ii);
                    
%                     disp(t2(ii-1))

                    
                    [yp, F, Q0] = fiber_dynamics_explicit_no_tendon_full(t2(ii-1),y(:,ii-1), parms);
                    
                    y(:,ii) = y(:,ii-1) + yp * dt;
                    
                    % for n, also need to shift
                    dx = yp(end-3) * dt;
                    xi = parms.xi - dx;
                    n = y(1:length(parms.xi),ii);
                    
                    nshift = interp1(parms.xi(:), n, xi(:));
                    nshift(isnan(nshift)) = 0;
                    nshift(nshift<0) = 0;
                    
                    y(1:length(parms.xi),ii) = nshift;
                    
                    t2(ii) = t2(ii-1) + dt;
                    
                    Qs = trapz(xi(:), [nshift(:) xi(:).*nshift(:)]);
                    F2(ii) = sum(Qs);
%                     F2(ii) = F + dx * Q0;

                end
                
                Ntot(id,j,f) = N;
                
            else
                
                % method of characteristics
                parms.method_of_characteristics = 1;
                
                Y = [];
                t1 = [];
   
                Nit = nan(1,(length(nzi)-1));
                for p = 1:(length(nzi)-1)
%                     disp(p)
                    
                    parms.vts = vs(nzi(p));
                    
                    sol = ode45(@(t,y,yp) fiber_dynamics_explicit_no_tendon_full(t,y, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, []);
                    
                    Y = [Y sol.y(:,1:end-1)];
                    t1 = [t1 sol.x(1:end-1)];
                    
                    X0 = sol.y(:,end);
                    
                    Nit(p) = sol.stats.nfevals;
                    
                end
                
                Ntot(id,j,f) = sum(Nit);
            end
            
            tsim(id,j,f) = toc;
            
            % compute force
            if f == 1
                Fi(:,id,j,f) = interp1(t2, F2, tlin);
                
            else
                
                L = Y(end-3,:);
                F1 = nan(1, size(Y,2));
                n = Y(1:end-4,:);

                for i = 1:size(Y,2)
                    xi = parms.xi + L(i);
                    F1(i) = trapz(xi, xi .* n(:,i)') + trapz(xi, n(:,i)');
                end

                Fi(:,id,j,f) = interp1(t1, F1, tlin);
            end
            
        end
    end
end

%% Gaussian-approximated
X0 = parms.x0';

Nit = nan(1,(length(nzi)-1));

Y = [];
t1 = [];
tic

for p = 1:(length(nzi)-1)

    parms.vts = vs(nzi(p));

    sol = ode45(@(t,y,yp) fiber_dynamics_explicit_no_tendon(t,y, parms), [aTs(nzi(p)) aTs(nzi(p+1))], X0, []);

    Y = [Y sol.y(:,1:end-1)];
    t1 = [t1 sol.x(1:end-1)];

    X0 = sol.y(:,end);

    Nit(p) = sol.stats.nfevals;

end

Napi = sum(Nit);
tapi = toc;

Fap = Y(1,:) + Y(2,:);

Fapi = interp1(t1, Fap, tlin);

%% forces
figure(2)

% subplot(221)
plot(tlin,Fi(:,end,1,1)); hold on
% plot(tlin,Fi(:,5,2,2),'--')

plot(tlin, Fapi,':')

% subplot(222)
% semilogy(tlin, abs(Fi - Fi(:,end)));

%%
close all
titles = {'Fixed strain range', 'Fixed bin width'};

figure(3)

color = get(gca,'colororder');

for j = 1:2
    for f = 1:2
      
        subplot(3,2,j)
        semilogx(Ns, tsim(:,j,f) ./ Ntot(:,j,f),'ko', 'markerfacecolor', color(f,:)); hold on
        box off
        xlabel('# bins')
        ylabel('Cost (s)')

        title(titles{j})
        subtitle('Cost')

        xlim([Ns(1) Ns(end-1)])
        
        ylim([0 1.2e-4])
%         
        if f == 2
            yline(tapi/Napi, 'k--', 'linewidth', 1)
        end
        
        subplot(3,2,j+2)
        loglog(Ns, mean(abs(Fi(:,:,j,f) - Fi(:,end,j,f)), 'omitnan'),'ko', 'markerfacecolor', color(f,:)); hold on
        box off
        subtitle('Error')
        xlabel('# bins')
        ylabel('Error (-)')

        xlim([Ns(1) Ns(end-1)])
        ylim([1e-10 1e0])
        
        if f == 2
            yline(mean(abs(Fapi(:) - Fi(:,end,j,f)), 'omitnan'),'k--', 'linewidth', 1)
        end
        
        subplot(3,2,j+4)
        semilogy(tsim(:,j,f) ./ Ntot(:,j,f), mean(abs(Fi(:,:,j,f) - Fi(:,end,j,f)), 'omitnan'),'ko', 'markerfacecolor', color(f,:)); hold on
        box off
        subtitle('Error versus cost')
        xlabel('Cost (s)')
        ylabel('Error (-)')
        
        if f == 2
            plot(tapi/Napi, mean(abs(Fapi(:) - Fi(:,end,j,f)), 'omitnan'), 'ko', 'markerfacecolor', [0 0 0])
            yline(mean(abs(Fapi(:) - Fi(:,end,j,f)), 'omitnan'),'k--', 'linewidth', 1)
            xline(tapi/Napi, 'k--', 'linewidth', 1)
        end
        
        ylim([1e-10 1e0])
        xlim([0 1.2e-4])
%         xlim([Ns(1) Ns(end-1)])

    end
end

%%

subplot(324)
legend('Tradit.', 'Charac.', 'Approx.', 'location', 'best')
legend boxon


subplot(326)
legend('Tradit.', 'Charac.', 'Approx.', 'location', 'best')
legend boxon
