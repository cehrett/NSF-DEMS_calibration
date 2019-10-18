%scratch


%% Try DCTO code for CTO example problem
clc ; clearvars -except dpath ; close all ;

% Set model observations
load([dpath,'Example\Ex_results\'...
'2019-10-16-raw_dat-3-6-6']);

sim_x = raw_dat.sim_xt(:,1);
sim_t = raw_dat.sim_xt(:,2:3);
sim_y = raw_dat.sim_y;
clear raw_dat;

% Settings
obs_var = 0.05; % observation error
obs_x_size  = 3;

% Set number of draws, burn_in for each mcmc:
M = 6e3; b = .5 ; burn_in = M*b;

% Define inputs mins and ranges 
xmin = 1.95;
xrange = .1;
tmin = [0 0];
trange = [3 6];

% Define comp. model & truth (version w/ nonstandardized output)
model_fn_ns = @(xt) Ex_sim(xt .* [xrange trange] + [xmin tmin]);

% Get mean and std using comp. model, define standardized version
X=lhsdesign(1000,3); 
Y=model_fn_ns(X);
mean_y = mean(Y) ; std_y = std(Y) ;
model_fn = @(xt) (model_fn_ns(xt) - mean_y)./std_y;

% Get initial design
obs_x = linspace(0,1,obs_x_size)' * xrange + xmin;
obs_y = 0 * ones(size(obs_x,1),1);


% Emulator mean
mean_sim = @(a,varargin) 0 ; 

% Get settings for DCTO
settings = MCMC_dual_calib_settings(sim_x,sim_t,[],sim_y,...
    obs_x,[],obs_y,...
    [],zeros(0,3),'min_x',xmin,'range_x',xrange,...
    'min_t1',tmin,'range_t1',trange,...
    'min_t2',[],'range_t2',[],...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',false,...
    'emulator_use',true,...
    'EmulatorMean',mean_sim,...
    'des_discrep',false,...
    'obs_var',obs_var,...
    'obs_discrep_use_MLEs',true,...
    'doplot',true);

% Perform dual calibration
% We need a loop because sometimes an ill-conditioned matrix early
% in the burn-in makes the whole calibration fail.
count = 0 ; err_count = 0 ; 
while count == err_count
    try
        res = MCMC_dual_calib(settings);
    catch ME
        warning('Warning: calibration failed. Retrying...');
        err_count = err_count + 1;
    end
    count = count + 1;
    if count >= 10, rethrow(ME) ; end
end