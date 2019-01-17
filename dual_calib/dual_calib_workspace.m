% Dual calibration workspace

%% Set path string and add paths
clc; clear all; close all;

direc = pwd; 
if direc(1)=='C' 
    dpath = 'C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc; 

% Add paths
addpath(dpath);
addpath([dpath,'stored_data']);
addpath([dpath,'dual_calib']);
addpath([dpath,'Example']);

%% Explore output of example function
clc ; clearvars -except dpath ; close all ;

% Define inputs
xmin = .5;
xrange = .5;
x = linspace(0,1);
t1min = 1.5;
t1range = 3;
t1=linspace(0,1);
t2min = 0;
t2range = 5;
t2 = linspace(0,1);

[X,T1,T2] = meshgrid(x,t1,t2) ; 
Y = reshape(dual_calib_example_fn(X(:),xmin,xrange,T1(:),t1min,t1range,...
    T2(:),t2min,t2range,0,1),length(x),length(t1),length(t2));

% Take a look
xidx=100;
xx=reshape(X(:,xidx,:),100,100);
tt1=reshape(T1(:,xidx,:),100,100);
tt2=reshape(T2(:,xidx,:),100,100);
surf(tt1 * t1range + t1min,tt2+t2range+t2min,reshape(Y(:,xidx,:),100,100));

%% Gather data to use for emulator and for "real" observations
clc ; clearvars -except dpath ; close all ;

% Define inputs mins and ranges 
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;

% Get simulation design points and observations
X = lhsdesign(250,3);
sim_x = X(:,1);
sim_t1 = X(:,2);
sim_t2 = X(:,3);
% Get simulation observations
sim_y = dual_calib_example_fn(...
    sim_x,xmin,xrange,sim_t1,t1min,t1range,sim_t2,t2min,t2range,0,1);


% Now set "real" theta1 and get "real" observations
theta1 = 2;
X = lhsdesign(30,2);
obs_x = X(:,1);
obs_t2 = X(:,2);
% Make a col vector based on true theta1
obs_t1 = (ones(size(obs_x,1),1) * theta1 - t1min)./t1range;
% Get "real" observations without noise or discrepancy
obs_y_noiseless = dual_calib_example_fn(...
    obs_x,xmin,xrange,obs_t1,t1min,t1range,obs_t2,t2min,t2range,0,1,false);

% Now noise it up (still keep no discrepancy for now)
sigma = sqrt(0.05); % This is noise s.d. of STANDARDIZED observations
std_true = std(sim_y);
obs_y = obs_y_noiseless + randn(size(obs_x,1),1) * sigma * std_true;


% Now set desired observations
des_x = linspace(0,1,8)';
des_y = zeros(size(des_x,1),1);


% Now package everything up and save it
% clear X std_true sigma obs_y_noiseless dpath
% save(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
%     'state-aware\state-aware_results\'...
%     '2019-01-15_dual_calib_raw_data']);


%% Get MLEs for covariance parameters for discrepancy GP
clc ; clearvars -except dpath ; close all ; 

% Load the raw data
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'state-aware\state-aware_results\'...
    '2019-01-15_dual_calib_raw_data']);

% unfinished


%% Perform dual calibration
% Here calibration will be performed without emulator or (true) discrepancy
clc ; clearvars -except dpath ; close all ;

% Load the raw data
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'state-aware\state-aware_results\'...
    '2019-01-15_dual_calib_raw_data']);
% Since we are not using emulator, empty out simulator observations
sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

% Get settings
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y);