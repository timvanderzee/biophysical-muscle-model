clear all; close all; clc
[username, githubfolder] = get_paths();

% th = [0 .05 .1 .25 .7 1.5];
th = [0 .07 .25 .7 1.5];
id = [3 1 7];

%%
figure(1)
visualize_data_SRS(th, id);

%%
% figure(1)
visualize_model_SRS_simple(th, id)

%%
set(gcf,'units','normalized','position',[.1 .1 .6 .3])