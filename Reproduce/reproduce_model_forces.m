clear all; close all; clc
cd(fileparts(which('reproduce_model_forces.m')));
cd ..
cd ..
githubfolder = cd;

% get folders
cd('biophysical-muscle-model')
addpath(genpath('Common'))
addpath(genpath('Reproduce'));
[~, modelfolder] = get_paths(githubfolder);
repofolder = cd;

%% step 0: create inputs (if they don't exist yet)
% the biophysical models have 3 inputs, namely:
% 1) calcium concentration
% 2) fiber length
% 3) fiber velocity
% note: the latter two depend on each other, because
% fiber velocity is the time-derivative of fiber length

cd(fullfile(githubfolder, 'biophysical-muscle-model', 'Reproduce', 'Process'))
visualize = 1;
make_ramp_protocols(githubfolder, modelfolder, visualize)

%% step 1: choose inputs
cd(fullfile(modelfolder, 'Protocols'))
[filename, pathname] = uigetfile('*.mat', 'Pick an input file');
load(fullfile(pathname, filename))

aTs = [0; Ts];
nzi = find(diff(aTs) > 0);

for i = 1:(length(nzi)-1)
    id = tis >= (aTs(nzi(i))) & tis <= (aTs(nzi(i+1)));
    
    input(i).t      = [aTs(nzi(i)) aTs(nzi(i+1))] - aTs(nzi(i));
    input(i).Cas    = mean(Cas(id)) * ones(1,2);
    input(i).L      = [Lis(find(id, 1, 'first')) Lis(find(id, 1, 'last'))];
    input(i).v      = vis(find(id, 1, 'last')) * ones(1,2);

end

% show length and velocity traces
if ishandle(1), close(1); end

color = lines(3);
for i = 1:length(input) % loop over phases
    t0 = aTs(nzi(i));
    
    figure(1)
    subplot(411)
    plot(input(i).t+t0, input(i).Cas, 'color', color(1,:), 'linewidth', 1.5); hold on; box off
    ylabel('[Ca]^{2+} (\muM)')
    xline(t0,'k--')
    title('Calcium concentration')
    ylim([0 35])
    
    subplot(412)
    plot(input(i).t+t0, input(i).L, 'color', color(1,:), 'linewidth', 1.5); hold on; box off
    ylabel('\DeltaFiber length (L_0)')
    xline(t0,'k--')
    title('Fiber length change')
    
    subplot(413)
    plot(input(i).t+t0, input(i).v, 'color', color(1,:), 'linewidth', 1.5); hold on; box off
    ylabel('Fiber velocity (L_0/s)')
    xline(t0,'k--')
    title('Fiber velocity')
    
end

figure(1)
subplot(411)
plot(tis, Cas, 'r:'); hold on

subplot(412)
plot(tis, Lis, 'r:'); hold on

subplot(413)
plot(tis, vis, 'r:'); hold on

%% step 2: specify model function
% you can choose between the following models:
% Hill-type SE, Hill-type no SE, 2-state XB, 2-state XB coop, 3-state XB coop, 4-state XB coop
model   = '2-state XB coop'; % see options above
odetype = 'explicit'; % type of differential equations
method  = 'approximated'; % solution method
% method = 'discretized'; % solution method
[modelfunc, odefunc, modelname] = look_up_model(model, odetype, method);

%% step 3: specify model parameters
% here we look up the parameters for a given fiber and model
cd(fullfile(repofolder, 'Reproduce', 'Parameters'))
fibername       = uigetdir('Pick a fiber');
parmsfilename    = fullfile(fibername, ['parms_', modelname, '.mat']);
load(parmsfilename, 'redparms')
newparms = complete_parms(redparms);

%% step 4: simulate model
[sol, out] = simulate_model(model, odefunc, modelfunc, input, newparms);

% visualize
figure(1)
subplot(414)

for i = 1:length(input)
    t0 = aTs(nzi(i));
    
    % plot the forces
    plot(out(i).t+t0, out(i).F, 'color', color(1,:), 'linewidth', 1.5); hold on; box off
    xline(t0,'k--')
end

xlabel('Time (s)')
ylabel('Force (F_0)')
title('Fiber force')

%% step 5: interpolate to input time and save
tsim = [];
Fsim = [];

for i = 1:length(input)
    t0 = aTs(nzi(i));
    
    tsim = [tsim out(i).t+t0];
    Fsim = [Fsim out(i).F];
end

% find the unique indices and interpolate
[tu, iu] = unique(tsim);
Fis = interp1(tu, Fsim(iu), tis);

figure(1)
subplot(414)
plot(tis, Fis, 'r:', 'linewidth', 1.5)
    
%% optional: loop over all conditions
% Hill-type SE, Hill-type no SE, 2-state XB, 2-state XB coop, 3-state XB coop, 4-state XB coop
model   = '3-state XB coop'; % see options above
odetype = 'explicit'; % type of differential equations
method  = 'approximated'; % solution method
[modelfunc, odefunc, modelname] = look_up_model(model, odetype, method);
visualize = 0;

if visualize
    if ishandle(1), close(1); end
    color = lines(1);

    figure(1)
    subplot(411)
    set(gca, 'Yscale', 'log')
    h(1) = line('xdata', [], 'ydata', [], 'color', color, 'linewidth', 1.5);
    xlim([0 14])
    ylim([1e-4 35])
    box off
    title('Calcium')

    subplot(412)
    h(2) = line('xdata', [], 'ydata', [], 'color', color, 'linewidth', 1.5);
    xlim([0 14])
    ylim([0 .08])
    box off
    title('Fiber length')

    subplot(413)
    h(3) = line('xdata', [], 'ydata', [], 'color', color, 'linewidth', 1.5);
    xlim([0 14])
    ylim([0 .5])
    box off              
    title('Fiber velocity')

    subplot(414)
    h(4) = line('xdata', [], 'ydata', [], 'color', color, 'linewidth', 1.5);
    xlim([0 14])
    ylim([0 2])
    box off              
    title('Fiber force')
    xlabel('Time (s)')
end

% loop over fibers
ffolders = dir(fullfile(repofolder, 'Reproduce', 'Parameters'));
for u = 1:length(ffolders)-2

    if ffolders(u+2).isdir
        fibername = ffolders(u+2).name;
        parmsfilename = fullfile(repofolder, 'Reproduce', 'Parameters', fibername, ['parms_', modelname, '.mat']);
        disp(['Loading: ', parmsfilename])
        load(parmsfilename, 'redparms')
        cd(fullfile(repofolder, 'Reproduce', 'Parameters'))
        newparms = complete_parms(redparms);
        set(gcf, 'name', fibername)

         % loop over pCas
        for j = 1:length(folders)-2
            folders = dir(fullfile(modelfolder, 'Protocols'));
            foldername = folders(j+2).name;
            files = dir(fullfile(folders(j+2).folder, foldername));
            
            % loop over protocols
            for k = 1:length(files)-2
                filename = files(k+2).name;
                fullname = fullfile(files(k+2).folder, filename);
%                 disp(['Loading: ', fullname])
                load(fullname)

                aTs = [0; Ts];
                nzi = find(diff(aTs) > 0);

                % loop over phases
                input = struct();
                for i = 1:(length(nzi)-1)
                    id = tis >= (aTs(nzi(i))) & tis <= (aTs(nzi(i+1)));

                    input(i).t      = [aTs(nzi(i)) aTs(nzi(i+1))] - aTs(nzi(i));
                    input(i).Cas    = mean(Cas(id)) * ones(1,2);
                    input(i).L      = [Lis(find(id, 1, 'first')) Lis(find(id, 1, 'last'))];
                    input(i).v      = vis(find(id, 1, 'last')) * ones(1,2);

                end

                % simulate model
                [sol, out] = simulate_model(model, odefunc, modelfunc, input, newparms);

                tsim = [];
                Fsim = [];
                for i = 1:length(input)
                    t0 = aTs(nzi(i));

                    tsim = [tsim out(i).t+t0];
                    Fsim = [Fsim out(i).F];
                end

                % find the unique indices and interpolate
                [tu, iu] = unique(tsim);
                Fis = interp1(tu, Fsim(iu), tis);
           
                if visualize
                    subplot(411)
                    set(h(1), 'xdata', tis, 'ydata', Cas)

                    subplot(412)
                    set(h(2), 'xdata', tis, 'ydata', Lis)

                    subplot(413)
                    set(h(3), 'xdata', tis, 'ydata', vis)

                    subplot(414)
                    set(h(4), 'xdata', tis, 'ydata', Fis)
                    drawnow
                end
                
                if contains(modelname, 'Hill')
                    subfolder = 'Hill';
                else
                    subfolder = method;
                end
            
                % save output
                savename = fullfile(modelfolder, subfolder, modelname, fibername, foldername, [fibername, '_', filename]);
                disp(['Saving: ', savename])
                save(savename, 'tis', 'Fis')
            end
        end
    end
end
