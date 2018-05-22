% Simulation example figures
% Here is code for figures based on the toy simulation example for the
% calibration to desired observations


clc; clear all; close all;

%% Set path string
direc = pwd; if direc(1)=='C' 
    dpath = 'C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc;

%% Add paths
addpath(dpath);
addpath([dpath,'stored_data']);
addpath([dpath,'Example']);
addpath([dpath,'Example\Ex_results']);

%% Get deflection, rotation "bands" as function of cost, with true trade-off
% First using the one big calibration including cost
load([dpath,'Example\Ex_results\'...
'2018-05-17_d0_incl_min_cost'],...
'results'); % This is the results
load([dpath,'Example\Ex_results\'...
    '2018-05-18_full_calib_emulator_output_estimates'],...
    'emout'); % This is the GP estimate of output at each draw of the MCMC


x=linspace(15,30); % cost values at which to plot

% Get quantiles of MCMC samples at each whole dollar amount
costs = 15:30;
data_at_cost = cell(length(costs),1);
medians = zeros(length(costs),2);
uqs = zeros(length(costs),2);
lqs = zeros(length(costs),2);
for ii = 1:16
    
    data = emout.output_means;
    data_at_cost{ii} = data(round(data(:,3))==costs(ii),1:2);
    medians(ii,:) = median(data_at_cost{ii});
    uqs(ii,:)     = quantile(data_at_cost{ii},1-alpha/2);
    lqs(ii,:)     = quantile(data_at_cost{ii},alpha/2);
    
end
plot(costs,medians(:,1))


% user defined values
alpha = 0.1;
k=1200;

% Get median of MCMC samples and plot
median_vals = nearquant(.5,x,emout.output_means(:,1),...
    emout.output_means(:,3),k);
plot(x,median_vals,'-o');

% Get median of direct data and plot