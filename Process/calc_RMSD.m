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

% mcodes = [1 1 1; 1 1 1; 1 1 3; 2 1 1];
mcodes = [1 1 1];

for kk = 1:size(mcodes,1)
    mcode = mcodes(kk,:);
    
    [output_mainfolder, filenames{kk}, opt_types{kk}, ~] = get_folder_and_model(mcode);
end

%% load data and parameters
fibers = {'12Dec2017a','13Dec2017a','13Dec2017b','14Dec2017a','14Dec2017b','18Dec2017a','18Dec2017b','19Dec2017a','6Aug2018a','6Aug2018b','7Aug2018a'};
iFs = [1 2 3, 5, 6, 7, 8, 10, 11];
pCas = repmat([4.5000    6.1000    6.2000    6.3000    6.4000    6.6000    9.0000],1,1,3);


% chosen ISIs, AMPs and pCas
ISIs = repmat([.001 .100 .316 1], 4,1,7);
AMPs = repmat([0 .0038 .0121 .0383]', 1, 4,7);

for i = 1:size(ISIs,3)
    ACTs(:,:,i) = pCas(i) * ones(size(ISIs,1), size(ISIs,2));
end

%%
% pCas = 6.1 * ones(size(AMPs));
RMSD = nan(size(ISIs,1), size(ISIs,2), length(iFs), size(mcodes,1));

for iF = 1:length(iFs)
    cd([output_mainfolder{2},'\data'])
    load([fibers{iFs(iF)},'_cor_new.mat'],'data');
    
    disp([fibers{iFs(iF)}])
    vs = {'\', '\'};
    
    for kk = 1:size(mcodes,1)
        disp(filenames{kk})
        output_folder = [opt_types{kk},'\normalized\with_PE_optimized\2_trials'];
        
        output_dir = [output_mainfolder{1}, '\', filenames{kk},vs{1}, output_folder];
        cd(output_dir)
        
        load([filenames{kk},'_F', num2str(iFs(iF)),'_best.mat'],'parms','exitflag','fopt','C0','Cbounds','model','P0','P')
        C = p_to_c(P, Cbounds);
        parms = C_to_parms(C, parms, parms.optvars);
        parms = calc_dependent_parms(parms);
        
        parms.act = 1;
        parms.Noverlap = 1;
        parms.vF_func =  @(vcerel,parms)parms.e(1)*log((parms.e(2)*vcerel./parms.vmax+parms.e(3))+sqrt((parms.e(2)*vcerel./parms.vmax+parms.e(3)).^2+1))+parms.e(4);
        
        Parms{kk} = parms;
        
        %     if kk == 1
        %         cd([mainfolder, '\GitHub\biophysical-muscle-model\Parameters'])
        %         load('parms_v5.mat')
        %         Parms{kk} = sparms(k);
        %     end
    end
    
    %% evaluate
    odeopt = odeset('maxstep', 1e-2);
    gamma = 108.3333; % length scaling
    
    for j = 1:size(ISIs,1)
        for i = 1:size(ISIs,2)        
            for k = 1:size(ISIs,3)
                
                % get data            
                [texp, Lexp, Fexp, Tsrel] = get_data(data, ISIs(j,i,k), AMPs(j,i,k), ACTs(j,i,k));
                
                % evaluate fit
                Ca = 10.^(-ACTs(j,i,k)+6);
                
                dTt = .0383/.4545; % test stretch (= constant)
                dTc = AMPs(j,i,k) / .4545; % conditioning stretch
                ISI = ISIs(j,i,k);
                tiso = 3;
                
                [tis, Cas, Lis, vis, ts] = create_input(tiso, dTt, dTc, ISI, Ca, 2000);
                
                % center around 2nd stretch
                t = tis + 3*dTt - tiso;
                
                
                for kk = 1:size(mcodes,1)
                    
                    parms = Parms{kk};
                    
                    if contains(filenames{kk}, 'Hill')
                        x0 = 0;
                    else
                        x0 = [parms.x0(2:end)'; 0];
                    end
                    
                    xp0 = zeros(size(x0));
                    
                    parms.ti = tis;
                    parms.vts = vis;
                    parms.Cas = Cas;
                    parms.Lts = Lis;
                    
                    % run simulation
                    if contains(filenames{kk}, 'Hill')
                        sol = ode15i(@(t,y,yp) hill_type_implicit(t,y,yp, parms), [0 max(tis)], x0, xp0, odeopt);
                    else
                        sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 max(tis)], x0, xp0, odeopt);
                    end
                    
                    
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
                    id = t < .16 & t > (-ISI - 2 * dTc - .1);
                    oFii = interp1(t(id), oFi(id), texp);
                    
                    RMSDs = sqrt((oFii - Fexp).^2) * 100;
                    
                    RMSD(j,i,iF,kk,k) = sqrt(mean((oFii - Fexp).^2, 'omitnan'));
                end
            end
        end
    end
end

return
%%

RMSDm = mean(RMSD, 5, 'omitnan');
RMSDi = nan(size(RMSDm));
RMSDc = nan(size(RMSDm));



% correct for individual offset
for iF = 1:9
    RMSDc(:,:,iF,:) = RMSDm(:,:,iF,:) - mean(RMSDm(:,:,iF,:),'all','omitnan') + mean(RMSDm(:),'all','omitnan');
end

return

%%

close all
figure(1)

% for i = 1:size(AMPs,3)
%     nexttile
i = 1;

surf(AMPs(:,:,i), ISIs(:,:,i), mean(RMSDc,3,'omitnan'));

set(gca,'YScale','log')
xlabel('Amplitude')
ylabel('Recovery time')

% title(['pCa = ', num2str(mean(ACTs(:,:,i),'all'))])

hold on

plot3(AMPs(1,1,i), ISIs(1,1,i), mean(RMSDc(1,1,:),3,'omitnan'),'ro')
plot3(AMPs(4,2,i), ISIs(4,2,i), mean(RMSDc(4,2,:),3,'omitnan'),'ro')

view(45, 45)
zlim([0 .05])
% end

%%

close all

ms = [4 4 2 2];
ns = [1 4 4 1];

titles = {'Overall','Short recovery','Long recovery','Long recovery','Short recovery'}; 
subtitles = {'','Large amplitude','Large amplitude','Small amplitude','Small amplitude'}; 

mnames = {'Hill','No coop.','Coop.', 'Coop. + FD'};

pcolors = flip(parula(7));

figure(1)
color = get(gca,'colororder');
colors = [color(2,:); pcolors(4:end-1,:)];

jid = 1;
jr = 1;

X = nan(length(jid), size(RMSDc,3));
Y = nan(length(jid), size(RMSDc,3));

mm = 1;

for i = 1:4
    subplot(1,5,i+1)
    
%     v = violinplot(squeeze(RMSDi(ms(i),ns(i),:,jid) - RMSDi(ms(i),ns(i),:,jr))*100, mnames, 'ViolinColor', colors,'Edgecolor', [1 1 1]);
    v = violinplot(squeeze(RMSDc(ms(i),ns(i),:,jid))*100, mnames, 'ViolinColor', colors,'Edgecolor', [1 1 1],'ShowMean',true, 'ShowMedian',false);
    
    for j = 1:length(jid)
        X(j,:) = v(j).ScatterPlot.XData;
        Y(j,:) = v(j).ScatterPlot.YData;
    end
    
    plot(X, Y, ':','color',[.5 .5 .5])
    
end

%%
subplot(151)
% v = violinplot((mmRMSDi(:,jid)-mmRMSDi(:,jr))*100, mnames, 'ViolinColor', colors,'Edgecolor', [1 1 1]); 
v = violinplot((mmRMSDc(:,jid))*100, mnames, 'ViolinColor', colors,'Edgecolor', [1 1 1],'ShowMean',true, 'ShowMedian',false); 

for j = 1:length(jid)
    X(j,:) = v(j).ScatterPlot.XData;
    Y(j,:) = v(j).ScatterPlot.YData;
end

plot(X, Y, ':','color',[.5 .5 .5])

drawnow  

%
figure(mm)

for i = 1:5
    subplot(1,5,i)
%     axis([.5 4.5 -100 100])
    axis([.5 4.5 -90 60])
%         axis([.5 3.5 -4 4])
        box off
%     axis tight
    title(titles{i})
    subtitle(subtitles{i},'Fontweight','bold')
    ylabel('\DeltaRMSD (%)')
    
    yline(0,'k-','linewidth',1.5)
end

subplot(151)
% a = annotation('textarrow',[.21 .21],[.7 .52],'String','Hill-type model','Fontsize',8,'HeadWidth',5,'HeadLength',5);

% make nice
for i = 2:5
    subplot(1,5,i)
    ylabel('')
    set(gca,'YTickLabel', '')
    
end
%
figure(mm)
set(gcf,'units','normalized','position',[.2 .2 .6 .3])

ABC = {'A','B','C','D','E'};
for i = 1:5
    subplot(1,5,i)
    ym = max(get(gca,'ylim'));
    text(0,ym*1.2, ABC{i},'fontsize',16)
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