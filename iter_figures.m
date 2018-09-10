%% Set path string and add paths
clc; clear all; close all;

direc = pwd; if direc(1)=='C' 
    dpath = 'C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc; 

% Add paths
addpath(dpath);
addpath([dpath,'stored_data']);
addpath([dpath,'Example']);
addpath([dpath,'Example\Ex_results']);

%% Heatmaps for posteriors in 3-part iter. calib. on toy sim. example
clc ; clearvars -except dpath ; close all ; 

% Load results
load([dpath,'Example\Ex_results\'...
    '2018-09-10_iterative_calibration'],...
    'all_results');

des_obs = all_results{1}.settings.desired_obs;

h1 = calib_heatmap(des_obs, all_results{1}.samples_os);
pos = h1.Children(4).Position; dist = pos(4); sc=1/7;
set(h1.Children(4),'Position', pos + [-dist*sc -dist*sc dist*sc dist*sc]);
title('Calibration iteration 1/3, \lambda_\delta = 1/256');
xlabel('\theta_1');ylabel('\theta_2');

h2 = calib_heatmap(des_obs, all_results{2}.samples_os);
pos = h2.Children(4).Position; dist = pos(4); sc=1/7;
set(h2.Children(4),'Position', pos + [-dist*sc -dist*sc dist*sc dist*sc]);
title('Calibration iteration 2/3, \lambda_\delta = 16.7');
xlabel('\theta_1');ylabel('\theta_2');

h3 = calib_heatmap(des_obs, all_results{3}.samples_os);
pos = h3.Children(4).Position; dist = pos(4); sc=1/7;
set(h3.Children(4),'Position', pos + [-dist*sc -dist*sc dist*sc dist*sc]);
title('Calibration iteration 3/3, \lambda_\delta = 1784.5');
xlabel('\theta_1');ylabel('\theta_2');

% Save them
set(h1,'Color','none');
set(h2,'Color','none');
set(h3,'Color','none');
export_fig('FIG_iter_calib_1o3','-png','-m3','-painters',h1);
export_fig('FIG_iter_calib_2o3','-png','-m3','-painters',h2);
export_fig('FIG_iter_calib_3o3','-png','-m3','-painters',h3);