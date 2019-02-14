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
surf(tt1 * t1range + t1min,tt2*t2range+t2min,reshape(Y(:,xidx,:),100,100));

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
sim_x = X(:,1) * xrange + xmin;
sim_t1 = X(:,2) * t1range + t1min;
sim_t2 = X(:,3) * t2range + t2min;
% Get simulation observations
sim_y = dual_calib_example_fn(sim_x,sim_t1,sim_t2);

% Now set "real" theta1 and get "real" observations
theta1 = 2;
X = lhsdesign(30,2);
obs_x = X(:,1) * xrange + xmin;
obs_t2 = X(:,2) * t2range + t2min;
% Make a col vector based on true theta1
obs_t1 = ones(size(obs_x,1),1) * theta1;
% Get "real" observations without noise or discrepancy
obs_y_noiseless = dual_calib_example_fn(obs_x,obs_t1,obs_t2);

% Now noise it up (still keep no discrepancy for now)
sigma = sqrt(0.05); % This is noise s.d. of STANDARDIZED observations
std_y = std(sim_y);
mean_y = mean(sim_y);
obs_y = obs_y_noiseless + randn(size(obs_x,1),1) * sigma * std_y;

% Now set desired observations
des_x = linspace(0,1,8)' * xrange + xmin;
des_y = zeros(size(des_x,1),1);

% Now package everything up and save it
clear X sigma obs_y_noiseless dpath
% save(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
%     'dual_calib\dual_calib_stored_data\'...
%     '2019-01-15_dual_calib_raw_data']);
% And get dpath back:
direc = pwd; 
if direc(1)=='C' 
    dpath = 'C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc; 


%% Get MLEs for covariance parameters for emulator GP
clc ; clearvars -except dpath ; close all ; 

% Load the raw data
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-01-15_dual_calib_raw_data']);

% Normalize the simulation inputs, standardize the outputs
x = (sim_x - xmin) ./ xrange ;
t1 = (sim_t1 - t1min) ./ t1range;
t2 = (sim_t2 - t2min) ./ t2range;
y = (sim_y - mean_y) ./ std_y;

% Define function for minimization
f = @(rl) ...
    -logmvnpdf(y',0,...
    gp_cov([rl(1) rl(2) rl(3)],[x t1 t2],[x t1 t2],rl(4),false) + ...
    1e-4*eye(size(x,1)));

% Perform minimization
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0 0 0 0];
ub = [1 1 1 Inf];
x0 = [.5 .5 .5 1];
[inp,fval,exitflag,output] = fmincon(f,x0,A,b,Aeq,beq,lb,ub);

% Result: rho values 0.9929, 0.7844, 0.0779; lambda value 0.0836.


%% Perform dual calibration without emulator or discrepancy
% Here calibration will be performed without emulator or (true) discrepancy
clc ; clearvars -except dpath ; close all ;

% Load the raw data
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-01-15_dual_calib_raw_data']);
% Since we are not using emulator, empty out simulator observations
sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

% Get settings
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',.5,'range_x',.5,...
    'min_t1',1.5,'range_t1',3,'min_t2',0,'range_t2',5,...
    'mean_y',mean_y,'std_y',std_y,'burn_in',.5,...
    'EmulatorCovHypers',[.5 .5 .5 Inf],'obs_discrep',false);

% Perform calibration
results = MCMC_dual_calib(settings);

% Save
% save(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
%     'dual_calib\dual_calib_stored_data\'...
%     '2019-01-25_dual_calib_nosim_nodiscr'],'results');

%% Perform dual calibration with emulator and without discrepancy
% Here calibration will be performed with emulator but without discrepancy
clc ; clearvars -except dpath ; close all ;

% Load the raw data
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-01-15_dual_calib_raw_data']);

% Get settings
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',.5,'range_x',.5,...
    'min_t1',1.5,'range_t1',3,'min_t2',0,'range_t2',5,...
    'mean_y',mean_y,'std_y',std_y,'burn_in',.5,...
    'obs_discrep',false);

% Perform calibration
results = MCMC_dual_calib(settings);

% Save
% save(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
%     'dual_calib\dual_calib_stored_data\'...
%     '2019-02-08_dual_calib_nodiscr'],'results');

%% Gather data for emulator and for "real" observations, with discrepancy
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
sim_x = X(:,1) * xrange + xmin;
sim_t1 = X(:,2) * t1range + t1min;
sim_t2 = X(:,3) * t2range + t2min;
% Get simulation observations
sim_y = dual_calib_example_fn(sim_x,sim_t1,sim_t2);

% Now set "real" theta1 and get "real" observations
theta1 = 2;
X = lhsdesign(30,2);
obs_x = X(:,1) * xrange + xmin;
obs_t2 = X(:,2) * t2range + t2min;
% Make a col vector based on true theta1
obs_t1 = ones(size(obs_x,1),1) * theta1;
% Get "real" observations without noise or discrepancy
obs_y_noiseless = dual_calib_example_fn(X(:,1),xmin,xrange,...
    (obs_t1-t1min)/t1range,t1min,t1range,X(:,2),t2min,t2range,0,1,true);

% Now noise it up (still keep no discrepancy for now)
sigma = sqrt(0.05); % This is noise s.d. of STANDARDIZED observations
std_y = std(sim_y);
mean_y = mean(sim_y);
obs_y = obs_y_noiseless + randn(size(obs_x,1),1) * sigma * std_y;

% Now set desired observations
des_x = linspace(0,1,8)' * xrange + xmin;
des_y = zeros(size(des_x,1),1);

% Now package everything up and save it
clear X sigma obs_y_noiseless dpath
% save(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
%     'dual_calib\dual_calib_stored_data\'...
%     '2019-02-14_dual_calib_raw_data']);
% And get dpath back:
direc = pwd; 
if direc(1)=='C' 
    dpath = 'C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc; 

%% Perform dual calibration without emulator but with discrepancy
% Here calibration will be performed without emulator or (true) discrepancy
clc ; clearvars -except dpath ; close all ;

% Load the raw data
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-02-14_dual_calib_raw_data']);
% Since we are not using emulator, empty out simulator observations
sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

% Get settings
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',.5,'range_x',.5,...
    'min_t1',1.5,'range_t1',3,'min_t2',0,'range_t2',5,...
    'mean_y',mean_y,'std_y',std_y,'burn_in',.5,...
    'EmulatorCovHypers',[.5 .5 .5 Inf],'obs_discrep',true);

% Perform calibration
results = MCMC_dual_calib(settings);

% Save
% save(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
%     'dual_calib\dual_calib_stored_data\'...
%     '2019-02-14_dual_calib_nosim_discr'],'results');

%% Perform dual calibration with emulator and discrepancy
% Here calibration will be performed with emulator but without discrepancy
clc ; clearvars -except dpath ; close all ;

% Load the raw data
load(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
    'dual_calib\dual_calib_stored_data\'...
    '2019-02-14_dual_calib_raw_data']);

% Get settings
settings = MCMC_dual_calib_settings(sim_x,sim_t1,sim_t2,sim_y,...
    obs_x,obs_t2,obs_y,des_x,des_y,'min_x',.5,'range_x',.5,...
    'min_t1',1.5,'range_t1',3,'min_t2',0,'range_t2',5,...
    'mean_y',mean_y,'std_y',std_y,'burn_in',.5,...
    'obs_discrep',true);

% Perform calibration
results = MCMC_dual_calib(settings);

% Save
% save(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\',...
%     'dual_calib\dual_calib_stored_data\'...
%     '2019-02-14_dual_calib'],'results');