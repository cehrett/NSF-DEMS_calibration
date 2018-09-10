% Iterative calibration with adaptive desired observations workspace

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

%% Perform 2-step calibration on toy sim example
clc ; clearvars -except dpath ; close all ;

% Load data
load([dpath,'Example\Ex_results\'...
'2018-05-28-raw_dat-3-12-12']);
sim_x = raw_dat.sim_xt(:,1);
sim_t = raw_dat.sim_xt(:,2:3);
sim_y = raw_dat.sim_y;
clear raw_dat;

% Get settings
desired_obs = [0 0 0 ] ; %[ 0.7130 0.7144 17.9220 ] ; 
settings = MCMC_settings(desired_obs,sim_x,sim_t,sim_y,...
    'Discrepancy',true,'M',2e4,'ObsVar','Constant',...
    'LambdaDeltaInit',1/16^2);

% Perform calibration
results = MCMC_discrepancy_true_fn(settings);
results.model_output.by_sample_true = ...
    Ex_sim([2*ones(size(results.samples,1),1) results.samples_os]);

all_results = cell(3,1);
all_results{1}=results;

% Get round 2 desired observation
[new_des_obs new_lambda_delta] = update_des_obs(results,desired_obs);

% Get round 2 settings
settings = MCMC_settings(new_des_obs,sim_x,sim_t,sim_y,...
    'Discrepancy',true,'M',2e4,'ObsVar','Constant',...
    'LambdaDeltaInit',new_lambda_delta);

% Perform calibration
results = MCMC_discrepancy_true_fn(settings);
results.model_output.by_sample_true = ...
    Ex_sim([2*ones(size(results.samples,1),1) results.samples_os]);

all_results{2}=results;

% Get round 3 desired observation
[new_des_obs new_lambda_delta] = update_des_obs(results,desired_obs);

% Get round 3 settings
settings = MCMC_settings(new_des_obs,sim_x,sim_t,sim_y,...
    'Discrepancy',true,'M',2e4,'ObsVar','Constant',...
    'LambdaDeltaInit',new_lambda_delta);

% Perform calibration
results = MCMC_discrepancy_true_fn(settings);
results.model_output.by_sample_true = ...
    Ex_sim([2*ones(size(results.samples,1),1) results.samples_os]);

all_results{3}=results;

% save results
% save([dpath,'Example\Ex_results\'...
%     '2018-09-10_iterative_calibration'],...
%     'all_results');

