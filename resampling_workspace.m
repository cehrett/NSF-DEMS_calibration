% Resampling workspace


clc; clear all; close all;

%% Add paths
addpath('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration');
addpath('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\stored_data');
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data');
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\');

load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\',...
    'Phase 1\stored_data\',...
    'results_3outputs_noVfThkPrior_1d',...
    '.mat']);

settings = MCMC_settings(size(results.samples,1)-1,...
    results.desired_obs,results.which_outputs);

[post_mean log_weights weights] = ...
    resample_mean(settings,results,[.7 .08 100])