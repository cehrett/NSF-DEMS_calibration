% Dynamic Vibration System DCTO analysis
% Here I analyze results obtained in dvs_dcto_workspace.m
% 2020-09-16
%% Set path string and add paths
clc; clear all; close all;

direc = pwd; if direc(1)=='C' 
    dpath = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc; 

% Add paths
addpath(dpath);
addpath([dpath,'stored_data']);
addpath([dpath,'dual_calib']);
addpath([dpath,'dual_calib\DVS_application']);
addpath([dpath,'dual_calib\DVS_application\data']);
addpath([dpath,'Example']);

% Change dir
cd(dpath);

%% Analyze inputs
clc ; clearvars -except dpath ; close all ; 

% Load the raw data
loadfile = [dpath, 'dual_calib\DVS_application\data',...
    '\2020-09-22_dvs_data_and_scaling_params'];
load(loadfile);

% Define true EM
obs_t1=6.2e10;

% Prepare figure with subplots
figure('units','normalized','outerposition',[0.5 0 .5 .35]);
num_ax = 3;
num_ax_cols = 3; %ceil(sqrt(num_ax)); % Number of columns
num_ax_rows = ceil(num_ax/num_ax_cols); % Number of rows
ax = gobjects(num_ax,1); % Preallocate axes array
for jj=1:num_ax
    ax(jj) = subplot(num_ax_rows,num_ax_cols,jj);
end

% Compare obs and sim x
subplot(ax(1));
histogram(sim_x,'Normalization','pdf','BinWidth',.05);
hold on;
histogram(obs_x,'Normalization','pdf','BinWidth',.05);
title('Mass');
ylim([0,12]);
legend('Simulated data','Experimental data');

% Compare obs and sim t1
subplot(ax(2));
histogram(sim_t1,'Normalization','pdf');
xline(obs_t1(1),'LineWidth',20,'color','red');
title('Elastic Modulus');

% Compare obs and sim t2
subplot(ax(3));
histogram(sim_t2,'Normalization','pdf','BinWidth',5);
hold on;
histogram(obs_t2,'Normalization','pdf','BinWidth',5);
title('Gain');
legend('Simulated data','Experimental data');

% Compare obs and sim y
% subplot(ax(4));
% histogram(sim_y,'Normalization','pdf','BinWidth',.05);
% hold on;
% histogram(obs_y,'Normalization','pdf','BinWidth',.05);
% title('Damping');

%% Produce model report using a set of DCTO results
clc ; clearvars -except dpath results ; close all ; 

% Load a specific file if no results in workspace
if ~exist('results')
    loadloc = [dpath,'dual_calib\DVS_application\data\',...
        '\2020-09-28_dvs_dcto_results'];
    load(loadloc);
end

% Define useful values
true_theta1 = 6.2e10;
x_min = results.settings.min_x;
x_range = results.settings.range_x;
t1_min = results.settings.min_t1;
t1_range = results.settings.range_t1;
t2_min = results.settings.min_t2;
t2_range = results.settings.range_t2;
y_mean = results.settings.mean_y;
y_std = results.settings.std_y;
sim_x_os = results.settings.sim_x * x_range + x_min;
sim_y_os = results.settings.sim_y * y_std + y_mean;
obs_x_os = results.settings.obs_x * x_range + x_min;
obs_y_os = results.settings.obs_y * y_std + y_mean;


% Gather some relevant parts of the results
burn_in = results.settings.burn_in;
theta1 = results.theta1((burn_in+1):end,:);
theta2 = results.theta2((burn_in+1):end,:);

% Prepare figure with subplots
figure('units','normalized','outerposition',[0.5 0 .5 1]);
num_ax = 4;
num_ax_cols = ceil(sqrt(num_ax)); % Number of columns
num_ax_rows = ceil(num_ax/num_ax_cols); % Number of rows
ax = gobjects(num_ax,1); % Preallocate axes array
for jj=1:num_ax
    ax(jj) = subplot(num_ax_rows,num_ax_cols,jj);
end

% Show histogram of calibration parameter with expected value
subplot(ax(1));
histogram(theta1);
xline(true_theta1,'r');
xlim([t1_min,t1_min+t1_range]);
title('Calibration input');

% Show histogram of design parameter
subplot(ax(2));
histogram(theta2);
xlim([t2_min,t2_min+t2_range]);
title('Design input');

%%% Show histogram of predicted model output along with sim, obs output
% First get model output
tic;
model_output_means = emulator_mean(results,(x_min+.5*x_range)*ones(size(theta1)),[theta1 theta2]);
toc;
% Now make histogram
subplot(ax(3));
histogram(model_output_means,'Normalization','pdf','BinWidth',.01);
hold on;
histogram(results.settings.sim_y * results.settings.std_y + results.settings.mean_y,'Normalization','pdf','BinWidth',.01);
histogram(results.settings.obs_y * results.settings.std_y + results.settings.mean_y,'Normalization','pdf','BinWidth',.01);
title('Model output');

%%% Show scatterplot of x values and model outputs
subplot(ax(4));
scatter(0.5 * ones(size(model_output_means)) * x_range + x_min, model_output_means);
hold on;
scatter(sim_x_os, sim_y_os);
scatter(obs_x_os, obs_y_os);
title('Output vs. x (Mass)');

%%% Print some useful info
fprintf('\n\n################################\n        MODEL REPORT\n################################\n');
% Print the GP mean
fprintf('GP Mean:');
disp(results.settings.mean_sim);

% Print whether obs_discrep used
fprintf('\bObservation discrepancy used:');
disp(results.settings.obs_discrep);

% Print whether modular
fprintf('\bModular:');
disp(results.settings.modular);

% Print calibration error
fprintf('\bCalibration error:');
disp(abs(mean(theta1)-true_theta1));

% Print Gelman-Rubin statistic
fprintf('\bGelman-Rubin statistic, theta1:');
disp(mcmcgr(reshape(results.theta1,1,1,size(results.theta1,1)),2));
fprintf('\bGelman-Rubin statistic, theta2:');
disp(mcmcgr(reshape(results.theta2,1,1,size(results.theta2,1)),2));

% Print posterior calibration input mean
fprintf('\bCalibration input mean:');
disp(mean(theta1));

% Print posterior design input mean
fprintf('\bDesign input mean:');
disp(mean(theta2));

% Print obs beta params for rho
fprintf('\bDesign input mean:');
disp(results.settings);

%% scr
size([results.theta1, results.theta2])


