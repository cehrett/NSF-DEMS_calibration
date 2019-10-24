%scratch


%% Try DCTO code for CTO example problem
clc ; clearvars -except dpath ; close all ;

% Set model observations
% load([dpath,'Example\Ex_results\'...
% '2019-10-16-raw_dat-3-6-6']);
% load([dpath,'Example\Ex_results\'...
% '2019-10-17-raw_dat-200obs']);
% load([dpath 'Example\Ex_results\'...
% '2019-10-21-raw_dat-30obs'],...
% 'raw_dat');
load([dpath,'Example\Ex_results\'...
'2019-10-21-raw_dat-300obs']);


sim_x = raw_dat.sim_xt(:,1);
sim_t = raw_dat.sim_xt(:,2:3);
sim_y = raw_dat.sim_y;
clear raw_dat;

% Settings
obs_var = 10000; % observation error var
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
% obs_y = repmat([0 0 0],size(obs_x,1),1); 
% obs_y = repmat([0.7130 0.7144 17.9220],size(obs_x,1),1);
obs_y = repmat(min(sim_y),obs_x_size,1);

% Emulator mean
mean_sim = @(a,varargin) repmat([0 0 0],size(a,1),1) ; 

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
    'obs_discrep_use_MLEs',false,...
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

%% Look at results on heatmap including true optimum location
% clc ; clearvars -except dpath ; close all;

% load true Pareto front
load([dpath,'Example\Ex_results\'...
    '2018-05-29_true_ctheta-output_nondom'],...
    'ctheta_output_nondom');

% Get standardized version
cntrl_mins = res.settings.min_x   ;
cntrl_rngs = res.settings.range_x ;
calib_mins = res.settings.min_t1   ;
calib_rngs = res.settings.range_t1 ;
output_mns = res.settings.mean_y       ;
output_sds = res.settings.std_y         ;

% Set desired observation
des_obs_os = [ 0 0 0 ] ; % on original scale
des_obs_os = obs_y(1,:) ; % on original scale
des_obs    = (des_obs_os - output_mns)./output_sds;

% Find closest point in Pareto front
true_pf_obs_os = ctheta_output_nondom(:,4:6); % Get outputs of true pf, os
true_pf_obs = (true_pf_obs_os - output_mns)./output_sds;
[m,idx] = min(sum((true_pf_obs - des_obs).^2,2));
optim = ctheta_output_nondom(idx,:);
optim_calib = optim(2:3);

% Get original scale samples
theta_os = res.theta1 ;

% Take a look against a heatmap
calib_heatmap(des_obs_os);
hold on;
burn_in = res.settings.burn_in;
scatter(theta_os(burn_in:end,1),...
    theta_os(burn_in:end,2));
plot(mean(theta_os(burn_in:end,1)),...
    mean(theta_os(burn_in:end,2)),'.m','MarkerSize',15);   
plot(optim_calib(1),optim_calib(2),'.g');

%% Get experimental design and simulator results to train emulator
clc ; clearvars -except dpath ; close all ;

% Set number of samples
n = 300;

% Get min,ranges
xmin = 1.95; xrange = 0.1;
t1min = 0; t1range = 3;
t2min = 0; t2range = 6;

% Get design
X = lhsdesign(n,3);
sim_x = X(:,1) * xrange + xmin;
sim_t1 = X(:,2) * t1range + t1min;
sim_t2 = X(:,3) * t2range + t2min;

% Get output
sim_xt = [sim_x sim_t1 sim_t2];
sim_y = Ex_sim(sim_xt);

% Package it
raw_dat = struct('sim_xt',sim_xt,'sim_y',sim_y);

plot(sim_x,sim_t2,'.');

save([dpath 'Example\Ex_results\'...
'2019-10-21-raw_dat-' int2str(n) 'obs'],...
'raw_dat');