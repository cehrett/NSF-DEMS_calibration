% Results gathered for use in the CTO draft paper.
% This file was not used for results presented in the draft that was
% submitted to Technometrics. I created this file in order to gather
% results when revising the paper after the Technometrics revision. The
% primary purpose of making this new file is that I intend to gather
% results without using a discrepancy function.

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

%% Get simulation observations
clc ; clearvars -except dpath ; close all ;

% Set number of samples
n = 100;

% Get min,ranges
xmin = 1.95; xrange = 0.1;
t1min = 0; t1range = 3;
t2min = 0; t2range = 6;

% Get design
X = lhsdesign(100,3);
sim_x = X(:,1) * xrange + xmin;
sim_t1 = X(:,2) * t1range + t1min;
sim_t2 = X(:,3) * t2range + t2min;

% Get output
sim_xt = [sim_x sim_t1 sim_t2];
sim_y = Ex_sim(sim_xt);

% Package it
raw_dat = struct('sim_xt',sim_xt,'sim_y',sim_y);

% plot(sim_x,sim_t2,'.');

% save([dpath 'Example\Ex_results\'...
% '2019-10-17-raw_dat-' int2str(n) 'obs'],...
% 'raw_dat');