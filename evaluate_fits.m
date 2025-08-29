clear all; close all; clc

mainfolder = 'C:\Users\u0167448\Documents\GitHub';
% addpath(genpath([mainfolder, '\muscle-thixotropy']))
% addpath(genpath([mainfolder, '\casadi-windows-matlabR2016a-v3.5.5']))

% Import casadi libraries
import casadi.*; 

%% specify data
load('active_trials.mat', 'Fm')
iFs = [1 2 3, 5, 6, 7, 8, 10, 11];
iFs = [5, 6, 7];

iF = 6; % fiber number (note: #4 and #9 lack pCa = 4.5)
Ks = find(Fm(:,iF) > 0.1); % only consider active trials
n = 3; % ISI number
m = 7; % AMP number
tiso = 3; % isometric time (s)

%% load data
[Data, tis, Lis, vis, Cas, ts] = get_data(iF,n,m,Ks,tiso);

% interpolate force
Fis = interp1(Data.t, Data.F, tis);

figure(1 + iF*10)
subplot(411)
plot(Data.t, Data.Ca,'r.'); hold on
plot(tis, Cas, 'b', 'linewidth',1); 
box off

subplot(412)
plot(Data.t, Data.v,'r.'); hold on
plot(tis, vis,'b',  'linewidth',1); 
box off

subplot(413)
plot(Data.t, Data.L,'r.'); hold on
plot(tis, Lis,'b',  'linewidth',1); 
box off

subplot(414)
plot(Data.t, Data.F,'r.');
box off


%% get parameters
mcode = [1 1 1];
vs = {'\', '\'};

cd([mainfolder, '\muscle-thixotropy\new_model\get_variable'])
[output_mainfolder, filename, opt_type, ~] = get_folder_and_model(mcode);

disp(filename)
output_folder = [opt_type,'\normalized\with_PE_optimized\2_trials'];

output_dir = [output_mainfolder{1}, '\', filename,vs{1}, output_folder];
cd(output_dir)

% get parallel parameters from that fiber
load([filename,'_F', num2str(iF),'_best.mat'],'parms','exitflag','fopt','C0','Cbounds','model','P0','P')
C = p_to_c(P, Cbounds);
oparms = parms;
oparms = C_to_parms(C, oparms, oparms.optvars);
kpe = oparms.kpe;
Fpe0 = oparms.Fpe0;
oparms.act = 1;
oparms.Noverlap = 1;
oparms.Fscale = 1.5;

% get other parameters from other fiber
ii = 6; % fiber from which parameters are obtained
load([filename,'_F', num2str(ii),'_best.mat'],'parms','exitflag','fopt','C0','Cbounds','model','P0','P')
C = p_to_c(P, Cbounds);
parms = C_to_parms(C, parms, parms.optvars);
parms = calc_dependent_parms(parms);  

parms.act = 1;
parms.Noverlap = 1;
parms.Fscale = 1.5;
parms.kpe = kpe;
parms.Fpe0 = Fpe0;