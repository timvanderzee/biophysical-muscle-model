clear all; close all; clc
[username, githubfolder] = get_paths();

%%
figure(1)
visualize_data_SRS();

%%
% figure(1)
visualize_model_SRS_simple()

%%
set(gcf,'units','normalized','position',[.1 .1 .6 .3])