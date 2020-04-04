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
fprintf('Ready.\n');

%% TSE Get mins, ranges, means, stds from objective fn
clc ; clearvars -except dpath ; close all ;

% Define inputs mins and ranges 
xmin = 1.95;
xrange = .1;
tmin = [0 0];
trange = [3 6];

minfn1 = @(x) Ex_sim(x) * [1 ; 0 ; 0] ; 
minfn2 = @(x) Ex_sim(x) * [0 ; 1 ; 0] ; 
minfn3 = @(x) Ex_sim(x) * [0 ; 0 ; 1] ; 
maxfn1 = @(x) -Ex_sim(x) * [1 ; 0 ; 0] ; 
maxfn2 = @(x) -Ex_sim(x) * [0 ; 1 ; 0] ; 
maxfn3 = @(x) -Ex_sim(x) * [0 ; 0 ; 1] ; 

[~,min1] = fmincon(minfn1,[2 2 2],[],[],[],[],[xmin tmin],...
    [xmin tmin]+[xrange trange]);
[~,min2] = fmincon(minfn2,[2 2 2],[],[],[],[],[xmin tmin],...
    [xmin tmin]+[xrange trange]);
[~,min3] = fmincon(minfn3,[2 2 2],[],[],[],[],[xmin tmin],...
    [xmin tmin]+[xrange trange]);

[~,max1] = fmincon(maxfn1,[2 2 2],[],[],[],[],[xmin tmin],...
    [xmin tmin]+[xrange trange]);
[~,max2] = fmincon(maxfn2,[2 2 2],[],[],[],[],[xmin tmin],...
    [xmin tmin]+[xrange trange]);
[~,max3] = fmincon(maxfn3,[2 2 2],[],[],[],[],[xmin tmin],...
    [xmin tmin]+[xrange trange]);

ex_fn_mins = [min1 min2 min3];
ex_fn_maxs = -[max1 max2 max3];
ex_fn_rngs = ex_fn_maxs - ex_fn_mins;

X=rand(1e6,3) .* [xrange trange] + [xmin tmin];
Y = Ex_sim(X);
mean_y = mean(Y);
std_y = std(Y);

% Get the mean and std of the distance from the mean
dists = sqrt(sum((Y-mean_y).^2,2));
mean_dist = mean(dists);
std_dist = std(dists);

% Get mean and s.d. of the distance from the mean on standardized scale
Y_std = (Y-mean_y)./std_y ;
dists_std = sqrt(sum(Y_std.^2,2));
mean_dist_std = mean(dists_std);
std_dist_std = std(dists_std);

%% TSE Gather PCTO results, estimating var, no emulator
clc ; clearvars -except dpath ; close all ;

% Settings
obs_var = 5000; % initial guess of observation error var
obs_x_size  = 3;
doplot = true;
verbose = true;

% Set number of draws, burn_in for each mcmc:
M = 2e4; b = .5 ;

% Define inputs mins and ranges 
xmin = 1.95;
xrange = .1;
tmin = [0 0];
trange = [3 6];

% Define comp. model & truth (version w/ nonstandardized output)
model_fn_ns = @(xt) Ex_sim(xt .* [xrange trange] + [xmin tmin]);

% Objective function mean and output were found via brute force (elsewhere)
mean_y = [0.9264    0.7925   20.9993] ; 
std_y = [0.0602    0.0899    3.1927] ;
model_fn = @(xt) (model_fn_ns(xt) - mean_y)./std_y;

% Get target outcomes
obs_x = linspace(0,1,obs_x_size)' * xrange + xmin;
% The following line sets the target to be the utopia point of the system,
% which was located using fmincon on each objective function
obs_y = repmat([0.7311 0.6675 15],obs_x_size,1); 
% Other possibilities not used here:
% obs_y = repmat([0 0 0],size(obs_x,1),1); 
% obs_y = repmat([0.7130 0.7144 17.9220],size(obs_x,1),1);

% Emulator mean: set to true function (standardized scale)
mean_sim = @(x,t1,t2) model_fn([x t1 t2]);
% If we were using an emulator, we would set:
% mean_sim = @(a,varargin) repmat([0 0 0],size(a,1),1) ; 

% Get settings for DCTO
settings = MCMC_dual_calib_settings([],[],[],[],...
    obs_x,[],obs_y,...
    [],zeros(0,3),'min_x',xmin,'range_x',xrange,...
    'min_t1',tmin,'range_t1',trange,...
    'min_t2',[],'range_t2',[],...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',false,...
    'emulator_use',false,...
    'EmulatorMean',mean_sim,...
    'des_discrep',false,...
    'obs_var',obs_var,...
    'obs_discrep_use_MLEs',false,...
    'doplot',doplot,...
    'verbose',verbose,...
    'obs_var_est',false,...
    'CTO',true);

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

% save([dpath,'Example\Ex_results\'...
%     '2019-10-31_CTO_grid3_12_12_varest'],...
%     'res');

%% TSE Gather results without PCTO, estimating var, no emulator
clc ; clearvars -except dpath ; close all ;

% Settings
obs_var = 0.5; % initial guess of observation error var
obs_x_size  = 3;
doplot = true;
verbose = true;

% Set number of draws, burn_in for each mcmc:
M = 2e4; b = .5 ;

% Define inputs mins and ranges 
xmin = 1.95;
xrange = .1;
tmin = [0 0];
trange = [3 6];

% Define comp. model & truth (version w/ nonstandardized output)
model_fn_ns = @(xt) Ex_sim(xt .* [xrange trange] + [xmin tmin]);

% Objective function mean and output were found via brute force (elsewhere)
mean_y = [0.9264    0.7925   20.9993] ; 
std_y = [0.0602    0.0899    3.1927] ;
model_fn = @(xt) (model_fn_ns(xt) - mean_y)./std_y;

% Get target outcomes
obs_x = linspace(0,1,obs_x_size)' * xrange + xmin;
% The following line sets the target to be the utopia point of the system,
% which was located using fmincon on each objective function
obs_y = repmat([0.7311 0.6675 15],obs_x_size,1); 
% Other possibilities not used here:
% obs_y = repmat([0 0 0],size(obs_x,1),1); 
% obs_y = repmat([0.7130 0.7144 17.9220],size(obs_x,1),1);

% Emulator mean: set to true function (standardized scale)
mean_sim = @(x,t1,t2) model_fn([x t1 t2]);
% If we were using an emulator, we would set:
% mean_sim = @(a,varargin) repmat([0 0 0],size(a,1),1) ; 

% Get settings for DCTO
settings = MCMC_dual_calib_settings([],[],[],[],...
    obs_x,[],obs_y,...
    [],zeros(0,3),'min_x',xmin,'range_x',xrange,...
    'min_t1',tmin,'range_t1',trange,...
    'min_t2',[],'range_t2',[],...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',false,...
    'emulator_use',false,...
    'EmulatorMean',mean_sim,...
    'des_discrep',false,...
    'obs_var',obs_var,...
    'obs_discrep_use_MLEs',false,...
    'doplot',doplot,...
    'verbose',verbose,...
    'obs_var_est',true,...
    'CTO',true);

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

locstr = [dpath,'Example\Ex_results\'...
    '2019-11-06_CTO_noemulator_varest_noPCTO'];
% save(locstr,'res');

%% TSE Identify nearby target outcome (without PCTO)
clc ; clearvars -except dpath ; close all ;

% Define inputs mins and ranges 
xmin = 1.95;
xrange = .1;
tmin = [0 0];
trange = [3 6];

% Objective function mean and output were found via brute force (elsewhere)
mean_y = [0.9264    0.7925   20.9993] ; 
std_y = [0.0602    0.0899    3.1927] ;

% The following line sets the target to be the utopia point of the system,
% which was located using fmincon on each objective function
target_original = [0.7311 0.6675 15]; 
target_original_std = (target_original - mean_y)./std_y;

% We will find the output that minimizes the distance to the utopia point,
% using fmincon. All outputs on standardized scale.
minfn = @(x) sum(((Ex_sim(x)-mean_y)./std_y-target_original_std).^2);
[X,~] = fmincon(minfn,[2 2 2],[],[],[],[],[xmin tmin],...
    [xmin tmin]+[xrange trange]);

% Now we can get the feasible target
target_feasible = Ex_sim(X);
target_feasible_std = (target_feasible - mean_y)./std_y;

% Move to a target 1 s.d. away from the feasible target, in the direction
% of the original target
direction_vector = target_feasible_std - target_original_std;
direction_vector_normd = direction_vector / norm(direction_vector);

% This s.d. of distances of model outputs from the mean was found elsewhere
% by brute force (on standardized outputs):
std_dists =  0.5332;

% New target (on standarized scale and original scale)
target_new_std = target_feasible_std - std_dists * direction_vector_normd ;
target_new = target_new_std .* std_y + mean_y;

% Also get the Euclidean distances of the standardized targets from the
% feasible point
target_original_distance = ...
    sqrt(sum((target_original_std - target_feasible_std).^2));
target_new_distance = ...
    sqrt(sum((target_new_std - target_feasible_std).^2));


%% TSE Gather results using PCTO-chosen target, estimating var, no emulator
% The target used here was actually found via brute force, not PCTO.
clc ; clearvars -except dpath ; close all ;

% Settings
obs_var = 0.5; % initial guess of observation error var
obs_x_size  = 3;
doplot = true;
verbose = true;

% Set number of draws, burn_in for each mcmc:
M = 2e4; b = .5 ;

% Define inputs mins and ranges 
xmin = 1.95;
xrange = .1;
tmin = [0 0];
trange = [3 6];

% Define comp. model & truth (version w/ nonstandardized output)
model_fn_ns = @(xt) Ex_sim(xt .* [xrange trange] + [xmin tmin]);

% Objective function mean and output were found via brute force (elsewhere)
mean_y = [0.9264    0.7925   20.9993] ; 
std_y = [0.0602    0.0899    3.1927] ;
model_fn = @(xt) (model_fn_ns(xt) - mean_y)./std_y;

% Get target outcomes
obs_x = linspace(0,1,obs_x_size)' * xrange + xmin;
% The following line sets the target to be one standard deviation away from
% the feasible region, on the line connecting the feasible region to the
% utopia point of the system (which was itself found via fmincon on each
% objective).
obs_y = repmat([0.7591    0.7573   18.6686],obs_x_size,1); 
% Other possibilities not used here:
% obs_y = repmat([0 0 0],size(obs_x,1),1); 
% obs_y = repmat([0.7130 0.7144 17.9220],size(obs_x,1),1);

% Emulator mean: set to true function (standardized scale)
mean_sim = @(x,t1,t2) model_fn([x t1 t2]);
% If we were using an emulator, we would set:
% mean_sim = @(a,varargin) repmat([0 0 0],size(a,1),1) ; 

% Get settings for DCTO
settings = MCMC_dual_calib_settings([],[],[],[],...
    obs_x,[],obs_y,...
    [],zeros(0,3),'min_x',xmin,'range_x',xrange,...
    'min_t1',tmin,'range_t1',trange,...
    'min_t2',[],'range_t2',[],...
    'mean_y',mean_y,'std_y',std_y,'M',M,'burn_in',b,...
    'obs_discrep',false,...
    'emulator_use',false,...
    'EmulatorMean',mean_sim,...
    'des_discrep',false,...
    'obs_var',obs_var,...
    'obs_discrep_use_MLEs',false,...
    'doplot',doplot,...
    'verbose',verbose,...
    'obs_var_est',true,...
    'CTO',true);

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

locstr = [dpath,'Example\Ex_results\'...
    '2019-11-04_CTO_noemulator_varest_afterPCTO'];
% save(locstr,'res');



%% TSE Get true optimum for utopia point, using fmincon
clc ; clearvars -except dpath ; close all ;

% Input the utopia point (found by brute force elsewhere)
des_obs = [0.7311 0.6675 15];

% Objective function mean and output were found via brute force (elsewhere)
mean_y = [0.9264    0.7925   20.9993] ; 
std_y = [0.0602    0.0899    3.1927] ;

% Define function to standardize model outputs
std_fn = @(x) (x-mean_y)./std_y;

% Define function which, given model inputs, gives square distance of model
% output at that point from the utopia point
dist_fn = @(x) sum((std_fn(Ex_sim([2 x])) - std_fn(des_obs)).^2,2);

[X,FVAL] = fmincon(dist_fn,[1.5 3],[],[],[],[],[0 0],[3 6])

%% TSE Get posterior means and standard deviations
clc ; clearvars -except dpath ; close all ;

% Load results
locstr = [dpath,'Example\Ex_results\'...
    '2019-11-04_CTO_noemulator_varest_afterPCTO'];
load(locstr);
res_PCTO = res;
locstr = [dpath,'Example\Ex_results\'...
    '2019-11-06_CTO_noemulator_varest_noPCTO'];
load(locstr);
res_utopia = res;
burn_in = res.settings.burn_in;

mean_t1_PCTO = mean(res_PCTO.theta1(burn_in:end,1));
mean_t2_PCTO = mean(res_PCTO.theta1(burn_in:end,2));
std_t1_PCTO = std(res_PCTO.theta1(burn_in:end,1));
std_t2_PCTO = std(res_PCTO.theta1(burn_in:end,2));
mean_t1_utopia = mean(res_utopia.theta1(burn_in:end,1));
mean_t2_utopia = mean(res_utopia.theta1(burn_in:end,2));
std_t1_utopia = std(res_utopia.theta1(burn_in:end,1));
std_t2_utopia = std(res_utopia.theta1(burn_in:end,2));

%% WTA perform PCTO 
clc ; clearvars -except dpath ; close all ; 

%%% Load raw data and desired observation
load([dpath,'stored_data\'...
    'raw_dat']);
sim_x = raw_dat(:,1);
sim_t = raw_dat(:,2:3);
sim_y = raw_dat(:,4:6); 
clear raw_dat

% Settings
obs_var = 5e7; % initial guess of observation error var
obs_x_size  = 3; % Number of target observations
doplot = true;
verbose = true;

% Set emulator hyperparameter MLEs
% These were found via fmincon.
EmulatorCovHypers = [...
    0.723919426870610   0.710378993535344   0.999999999999718;
    0.978844923780470   0.972316694561107   0.998830926614838;
    0.990609430210837   0.988186044400943   0.998556385291986;
    0.017707123599447   0.026113351305598   0.000949737279367];

% Set number of draws, burn_in for each mcmc:
M = 1e4; b = .5 ; burn_in = M*b;

% Define inputs mins and ranges 
xmin = min(sim_x);
xrange = range(sim_x);
tmin = min(sim_t);
trange = range(sim_t);

% Get output means and stds
mean_y = mean(sim_y) ; std_y = std(sim_y) ;

% Get target design
obs_x = linspace(0,1,obs_x_size)' * xrange + xmin;
% Select the target design to be the utopia point, estimated from the
% existing observations of the computer model output
obs_y = repmat(min(sim_y),obs_x_size,1); % Utopia point


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
    'EmulatorCovHypers',EmulatorCovHypers,...
    'des_discrep',false,...
    'obs_var',obs_var,...
    'obs_discrep_use_MLEs',false,...
    'doplot',doplot,...
    'verbose',verbose,...
    'obs_var_est',false,...
    'CTO',true);

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

% save([dpath,'stored_data\'...
%     '2019-11-04_PCTO'],...
%     'res');

%% WTA Given PCTO results, get predicted output at sampled points
clc ; clearvars -except dpath ; close all ;

% Load the results
locstr = [dpath,'stored_data\'...
    '2019-11-04_PCTO'];
load(locstr);

% Define inputs mins and ranges 
xmin = res.settings.min_x;
xrange = res.settings.range_x;
tmin = res.settings.min_t1;
trange = res.settings.range_t1;

% Get posterior samples
burn_in = res.settings.burn_in;
t1 = res.theta1(:,1); t2 = res.theta1(:,2);

% Get set up for loop that will find posterior predictive dist
pred_xt_all = ([.5*ones(size(t1,1),1) t1 t2] - [0 tmin])./[1 trange];
sim_xt = [res.settings.sim_x res.settings.sim_t1 res.settings.sim_t2];
prior_mean_xt_all = res.settings.mean_sim(pred_xt_all) ; 
eta = res.settings.sim_y;
% num_mod_obs = 8
% mod_obs_idx = randsample(1:504,num_mod_obs);
% sim_xt = sim_xt(mod_obs_idx,:);
% eta = eta(mod_obs_idx,:);
prior_mean_simxt = res.settings.mean_sim(sim_xt) ; 

% Get posterior predictive distribution for each output
% Only do maxsimul samples at once, to ease computational expense
maxsimul = 1000;
for jj = 1 : ceil(size(t1,1)/maxsimul)
    startpt = (jj-1) * maxsimul + 1 ;
    endpt = min(jj * maxsimul, size(t1,1));
    pred_xt = pred_xt_all(startpt:endpt,:);
    for ii=1:3
        cov_pxt_xt = ...
            gp_cov(res.settings.emulator_rho(:,ii), pred_xt, sim_xt, ...
            res.settings.emulator_lambda(ii),false);
        cov_xt_xt = ...
            gp_cov(res.settings.emulator_rho(:,ii), sim_xt, sim_xt, ...
            res.settings.emulator_lambda(ii),false);
        cov_xt_xt = cov_xt_xt + 1e-5*eye(size(cov_xt_xt));
        cov_pxt_pxt = ...
            gp_cov(res.settings.emulator_rho(:,ii), pred_xt, pred_xt,...
            res.settings.emulator_lambda(ii),false);
        pred_mean_xt_current(:,ii) = ...
            prior_mean_xt_all(startpt:endpt,ii) + ...
            cov_pxt_xt * (cov_xt_xt \ ...
            (eta(:,ii) - prior_mean_simxt(:,ii))) ; 
        pred_cov_xt_current(:,:,ii) = ...
            cov_pxt_pxt - cov_pxt_xt * (cov_xt_xt \ cov_pxt_xt') ; 
        fprintf('\n%g of %g, %g of 3 Done\n',jj,...
            ceil(size(t1,1)/maxsimul),ii);
        pred_mean_xt(startpt:endpt,ii) = pred_mean_xt_current(:,ii);
        pred_cov_xt(startpt:endpt,startpt:endpt,ii) = ...
            pred_cov_xt_current(:,:,ii);
    end
end
% Put back on original scale
mean_y = res.settings.mean_y; std_y = res.settings.std_y;
pred_mean_xt_os = pred_mean_xt .* std_y + mean_y;
mean(pred_mean_xt_os)

% Take a look
whichidx=3;
histogram(pred_mean_xt_os(:,whichidx),'Normalization','pdf');

% Now resample from the relevant gaussian dists, to get full uncertainty
for ii = 1:3
    pred_stds_xt(:,ii) = sqrt(diag(pred_cov_xt(:,:,ii)));
    pred_samps_xt(:,ii) = randn(size(pred_mean_xt,1),1) .* ...
        pred_stds_xt(:,ii) + pred_mean_xt(:,ii);
end
% Put back on original scale
pred_samps_xt_os = pred_samps_xt .* std_y + mean_y;
mean(pred_samps_xt_os)

% Take a look
hold on;
histogram(pred_samps_xt_os(:,whichidx),'Normalization','pdf');


%%%%%%% Get prior predictive output 

% Number of samples
m = 1e4 ; 

% Get theta values
thetavals = rand(m,2) .* trange + tmin ; 
t1 = thetavals(:,1) ; t2 = thetavals(:,2); 

% Get set up for loop that will find prior predictive dist
pred_xt = ([.5*ones(size(t1,1),1) t1 t2] - [0 tmin])./[1 trange];
prior_mean_xt = res.settings.mean_sim(pred_xt) ; 
clear pred_mean_xt pred_cov_xt pred_stds_xt pred_samps_xt

% Get posterior predictive distribution for each output
pred_mean_xt = nan(m,3);
pred_cov_xt = nan(m,m,3);
for ii=1:3
    cov_pxt_xt = ...
        gp_cov(res.settings.emulator_rho(:,ii), pred_xt, sim_xt, ...
        res.settings.emulator_lambda(ii),false);
    cov_xt_xt = gp_cov(res.settings.emulator_rho(:,ii), sim_xt, sim_xt, ...
        res.settings.emulator_lambda(ii),false);
    cov_xt_xt = cov_xt_xt + 1e-5*eye(size(cov_xt_xt));
    % cov_xt_xt_inv = inv( cov_xt_xt ) ;  
    cov_pxt_pxt = ...
        gp_cov(res.settings.emulator_rho(:,ii), pred_xt, pred_xt,...
        res.settings.emulator_lambda(ii),false);
    pred_mean_xt(:,ii) = ...
        prior_mean_xt(:,ii) + cov_pxt_xt * (cov_xt_xt \ ...
        (eta(:,ii) - prior_mean_simxt(:,ii))) ; 
    pred_cov_xt(:,:,ii) = ...
        cov_pxt_pxt - cov_pxt_xt * (cov_xt_xt \ cov_pxt_xt') ; 
    fprintf('\n%g/3 done\n',ii);
end
% Put back on original scale
mean_y = res.settings.mean_y; std_y = res.settings.std_y;
prior_mean_xt_os = pred_mean_xt .* std_y + mean_y;
mean(prior_mean_xt_os)

% Take a look
histogram(prior_mean_xt_os(:,1),'Normalization','pdf');

% Now resample from the relevant gaussian dists, to get full uncertainty
for ii = 1:3
    pred_stds_xt(:,ii) = sqrt(diag(pred_cov_xt(:,:,ii)));
    pred_samps_xt(:,ii) = randn(size(pred_mean_xt,1),1) .* ...
        pred_stds_xt(:,ii) + pred_mean_xt(:,ii);
end
% Put back on original scale
prior_samps_xt_os = pred_samps_xt .* std_y + mean_y;
mean(prior_samps_xt_os)

% Take a look
hold on;
histogram(prior_samps_xt_os(:,1),'Normalization','pdf');

%%%%%%%%%%% Make a figure
% xlims = [0.7135    1.0765 ; 0.6425    1.0275 ; 14.2500   30.7500];
% ylims = [0   95 ; 0   120 ; 0    3.5000];
fig1 = figure('Position',[10 10 900 300],'Color','white');
for ii = 1:3
    subplot(1,3,ii);
    histogram(prior_mean_xt_os(:,ii),'Normalization','pdf',...
        'EdgeColor','none');
    hold on;
    histogram(pred_mean_xt_os(:,ii),'Normalization','pdf',...
        'EdgeColor','none');
    xlabel(sprintf('y_%g',ii));
    set(gca,'Yticklabel',[]);
    set(gca,'Ytick',[]);
%     xlim(xlims(ii,:));
%     ylim(ylims(ii,:));
end
% export_fig 'posterior_predictive_dists_without_GP_uncertainty' -m2

fig2 = figure('Position',[10 10 900 300],'Color','white');
for ii = 1:3
    subplot(1,3,ii);
    histogram(prior_samps_xt_os(:,ii),'Normalization','pdf',...
        'EdgeColor','none');
    hold on;
    histogram(pred_samps_xt_os(:,ii),'Normalization','pdf',...
        'EdgeColor','none');
    xlabel(sprintf('y_%g',ii));
    set(gca,'Yticklabel',[]);
    set(gca,'Ytick',[]);
%     xlim(xlims(ii,:));
%     ylim(ylims(ii,:));
end
% export_fig 'posterior_predictive_dists_with_GP_uncertainty' -m2

% Save the posterior predictive output
model_output.means = pred_mean_xt;
model_output.vars = cell2mat(...
    arrayfun(@(b)diag(pred_cov_xt(:,:,b)),1:3,'UniformOutput',false));
res.model_output = model_output;
% save(locstr,'res');

%% WTA Use PCTO results to find new target
clc ; clearvars -except dpath ; close all ;

% Load results from PCTO
locstr = [dpath,'stored_data\'...
    '2019-11-04_PCTO'];
load(locstr);

% Filter the results to keep only non-dominated ones
[PF, PFidx] = nondominated(res.model_output.means) ; 

% Convert back to original scale
mean_y = res.settings.mean_y ; std_y = res.settings.std_y;
PF_os = PF .* std_y + mean_y;

% Take a look
scatter3(PF(:,1),PF(:,2),PF(:,3))   ;
% Visually, we identify the point elbow as the
% "elbow" of the PF, which we select as a target. But we need to back away
% from the PF a bit: [0.76 0.091 130]
elbow_os = [0.7498 0.0902 140.3302];
elbow = (elbow_os - mean_y)./std_y;
% Previously: [0.74765, 0.089937, 142.4953], yielding [0.74 0.089 125]
% Previously: [0.74855, 0.090066, 140.9148], yielding [0.745 0.090 125]
% Direction to go from [0.7498 0.0902 140.3302]: 
% [-0.000800000000000   0.000100000000000 -5.330199999999991]
new_target = elbow - [.2 .2 1.1275];
new_target_os = new_target .* std_y + mean_y;
hold on;
scatter3(new_target(1),new_target(2),new_target(3))

% Verify that the closest point on the PF is the elbow we want:
dists = sum((PF - new_target).^2,2);
[~,idx]=min(dists);
near_point_on_PF = PF(idx,:);
scatter3(near_point_on_PF(1),near_point_on_PF(2),near_point_on_PF(3))
disp(near_point_on_PF-elbow)
% Yes, it targets the point we wanted.

%% WTA Perform CTO using target chosen in PCTO
clc ; clearvars -except dpath ; close all ; 

%%% Load raw data
load([dpath,'stored_data\'...
    'raw_dat']);
sim_x = raw_dat(:,1);
sim_t = raw_dat(:,2:3);
sim_y = raw_dat(:,4:6); 
clear raw_dat

% Settings
obs_var = 0.5; % initial guess of observation error var
obs_x_size  = 3; % Number of target observations
doplot = true;
verbose = true;

% Set emulator hyperparameter MLEs
% These were found via fmincon.
EmulatorCovHypers = [...
    0.723919426870610   0.710378993535344   0.999999999999718;
    0.978844923780470   0.972316694561107   0.998830926614838;
    0.990609430210837   0.988186044400943   0.998556385291986;
    0.017707123599447   0.026113351305598   0.000949737279367];

% Set number of draws, burn_in for each mcmc:
M = 1e4; b = .5 ;

% Define inputs mins and ranges 
xmin = min(sim_x);
xrange = range(sim_x);
tmin = min(sim_t);
trange = range(sim_t);

% Get output means and stds
mean_y = mean(sim_y) ; std_y = std(sim_y) ;

% Get target design
obs_x = linspace(0,1,obs_x_size)' * xrange + xmin;
% Select the target design to be the point chosen in PCTO (elsewhere)
obs_y = repmat(...
    [0.742497187630386   0.089184738492403  71.112749103024015],...
    obs_x_size,1);


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
    'EmulatorCovHypers',EmulatorCovHypers,...
    'des_discrep',false,...
    'obs_var',obs_var,...
    'obs_discrep_use_MLEs',false,...
    'doplot',doplot,...
    'verbose',verbose,...
    'obs_var_est',true,...
    'CTO',true);

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

locstr = [dpath,'stored_data\'...
    '2019-11-05_CTO'];
% save(locstr,'res');

%% WTA Get posterior predictive distribution from CTO after PCTO, and prior
clc ; clearvars -except dpath ; close all ;

% Load the results
locstr = [dpath,'stored_data\'...
    '2019-11-05_CTO'];
load(locstr);

% Define inputs mins and ranges 
xmin = res.settings.min_x;
xrange = res.settings.range_x;
tmin = res.settings.min_t1;
trange = res.settings.range_t1;

% Get posterior samples
burn_in = res.settings.burn_in;
t1 = res.theta1(:,1); t2 = res.theta1(:,2);

% Get set up for loop that will find posterior predictive dist
pred_xt_all = ([.5*ones(size(t1,1),1) t1 t2] - [0 tmin])./[1 trange];
sim_xt = [res.settings.sim_x res.settings.sim_t1 res.settings.sim_t2];
prior_mean_xt_all = res.settings.mean_sim(pred_xt_all) ; 
eta = res.settings.sim_y;
% num_mod_obs = 8
% mod_obs_idx = randsample(1:504,num_mod_obs);
% sim_xt = sim_xt(mod_obs_idx,:);
% eta = eta(mod_obs_idx,:);
prior_mean_simxt = res.settings.mean_sim(sim_xt) ; 
cov_xt_xt = nan(size(sim_xt,1),size(sim_xt,1),3);
for ii=1:3
    cov_xt_xt(:,:,ii) = ...
        gp_cov(res.settings.emulator_rho(:,ii), sim_xt, sim_xt, ...
        res.settings.emulator_lambda(ii),false);
    cov_xt_xt(:,:,ii) = cov_xt_xt(:,:,ii) + ...
        1e-5*eye(size(cov_xt_xt(:,:,ii)));
end

% Get posterior predictive distribution for each output
% Only do maxsimul samples at once, to ease computational expense
maxsimul = 1000;
for jj = 1 : ceil(size(t1,1)/maxsimul)
    startpt = (jj-1) * maxsimul + 1 ;
    endpt = min(jj * maxsimul, size(t1,1));
    pred_xt = pred_xt_all(startpt:endpt,:);
    for ii=1:3
        cov_pxt_xt = ...
            gp_cov(res.settings.emulator_rho(:,ii), pred_xt, sim_xt, ...
            res.settings.emulator_lambda(ii),false);
        cov_pxt_pxt = ...
            gp_cov(res.settings.emulator_rho(:,ii), pred_xt, pred_xt,...
            res.settings.emulator_lambda(ii),false);
        pred_mean_xt_current(:,ii) = ...
            prior_mean_xt_all(startpt:endpt,ii) + ...
            cov_pxt_xt * (cov_xt_xt(:,:,ii) \ ...
            (eta(:,ii) - prior_mean_simxt(:,ii))) ; 
        pred_cov_xt_current(:,:,ii) = ...
            cov_pxt_pxt - cov_pxt_xt * (cov_xt_xt(:,:,ii) \ cov_pxt_xt') ; 
        fprintf('\n%g of %g, %g of 3 Done\n',jj,...
            ceil(size(t1,1)/maxsimul),ii);
        pred_mean_xt(startpt:endpt,ii) = pred_mean_xt_current(:,ii);
        pred_cov_xt(startpt:endpt,startpt:endpt,ii) = ...
            pred_cov_xt_current(:,:,ii);
    end
end
% Put back on original scale
mean_y = res.settings.mean_y; std_y = res.settings.std_y;
pred_mean_xt_os = pred_mean_xt .* std_y + mean_y;
mean(pred_mean_xt_os)

% Now resample from the relevant gaussian dists, to get full uncertainty
for ii = 1:3
    pred_stds_xt(:,ii) = sqrt(diag(pred_cov_xt(:,:,ii)));
    pred_samps_xt(:,ii) = randn(size(pred_mean_xt,1),1) .* ...
        pred_stds_xt(:,ii) + pred_mean_xt(:,ii);
end
% Put back on original scale
pred_samps_xt_os = pred_samps_xt .* std_y + mean_y;
mean(pred_samps_xt_os)


%%%%%%% Get prior predictive output 

% Number of samples
m = 1e4 ; 

% Get theta values
thetavals = rand(m,2) .* trange + tmin ; 
t1 = thetavals(:,1) ; t2 = thetavals(:,2); 

% Get set up for loop that will find prior predictive dist
pred_xt = ([.5*ones(size(t1,1),1) t1 t2] - [0 tmin])./[1 trange];
prior_mean_xt = res.settings.mean_sim(pred_xt) ; 

% Get posterior predictive distribution for each output
pred_mean_xt_prior = nan(m,3);
pred_cov_xt_prior = nan(m,m,3);
for ii=1:3
    cov_pxt_xt = ...
        gp_cov(res.settings.emulator_rho(:,ii), pred_xt, sim_xt, ...
        res.settings.emulator_lambda(ii),false);
    % cov_xt_xt_inv = inv( cov_xt_xt ) ;  
    cov_pxt_pxt = ...
        gp_cov(res.settings.emulator_rho(:,ii), pred_xt, pred_xt,...
        res.settings.emulator_lambda(ii),false);
    pred_mean_xt_prior(:,ii) = ...
        prior_mean_xt(:,ii) + cov_pxt_xt * (cov_xt_xt(:,:,ii) \ ...
        (eta(:,ii) - prior_mean_simxt(:,ii))) ; 
    pred_cov_xt_prior(:,:,ii) = ...
        cov_pxt_pxt - cov_pxt_xt * (cov_xt_xt(:,:,ii) \ cov_pxt_xt') ; 
    fprintf('\n%g/3 done\n',ii);
end
% Put back on original scale
mean_y = res.settings.mean_y; std_y = res.settings.std_y;
prior_mean_xt_os = pred_mean_xt_prior .* std_y + mean_y;
mean(prior_mean_xt_os)


% Now resample from the relevant gaussian dists, to get full uncertainty
for ii = 1:3
    pred_stds_xt(:,ii) = sqrt(diag(pred_cov_xt_prior(:,:,ii)));
    pred_samps_xt(:,ii) = randn(size(pred_mean_xt_prior,1),1) .* ...
        pred_stds_xt(:,ii) + pred_mean_xt_prior(:,ii);
end
% Put back on original scale
prior_samps_xt_os = pred_samps_xt .* std_y + mean_y;
mean(prior_samps_xt_os)


%%%%%%%%%%% Make a figure
% xlims = [0.7135    1.0765 ; 0.6425    1.0275 ; 14.2500   30.7500];
% ylims = [0   95 ; 0   120 ; 0    3.5000];
fig1 = figure('Position',[10 10 900 300],'Color','white');
for ii = 1:3
    subplot(1,3,ii);
    histogram(prior_mean_xt_os(5000:end,ii),'Normalization','pdf',...
        'EdgeColor','none');
    hold on;
    histogram(pred_samps_xt_os(5000:end,ii),'Normalization','pdf',...
        'EdgeColor','none');
    xlabel(sprintf('y_%g',ii));
    set(gca,'Yticklabel',[]);
    set(gca,'Ytick',[]);
%     xlim(xlims(ii,:));
%     ylim(ylims(ii,:));
end
% export_fig 'posterior_predictive_dists_without_GP_uncertainty' -m2

% Save the posterior predictive output
model_output.means = pred_mean_xt;
model_output.vars = cell2mat(...
    arrayfun(@(b)diag(pred_cov_xt(:,:,b)),1:3,'UniformOutput',false));
res.model_output = model_output;
% save(locstr,'res');

prior_model_output.theta = thetavals;
prior_model_output.means = pred_mean_xt_prior;
prior_model_output.vars = cell2mat(...
    arrayfun(...
        @(b)diag(pred_cov_xt_prior(:,:,b)),1:3,'UniformOutput',false));
locstr2 = [dpath,'stored_data\'...
    '2019-11-06_prior_predictive_distributions'];
% save(locstr2,'prior_model_output')

%% WTA Perform CTO on bivariate objective over cost grid (for PF est)
clc ; clearvars -except dpath ; close all ;

%%% Load raw data
load([dpath,'stored_data\'...
    'raw_dat']);
sim_x = raw_dat(:,1);
sim_t = raw_dat(:,2:3);
sim_y = raw_dat(:,[4,6]); 
clear raw_dat

% Define inputs mins and ranges 
xmin = min(sim_x);
xrange = range(sim_x);
tmin = min(sim_t);
trange = range(sim_t);

% Settings
obs_var = 0.01; % initial guess of observation error var
obs_x_size  = 3; % Number of target observations
doplot = true;
verbose = true;

% Set emulator hyperparameter MLEs
% These were found via fmincon.
EmulatorCovHypers = [...
    0.723919426870610      0.999999999999718;
    0.978844923780470      0.998830926614838;
    0.990609430210837      0.998556385291986;
    0.017707123599447      0.000949737279367];

% Set number of draws, burn_in for each mcmc:
M = 1e4; b = .5 ;

% Get output means and stds
mean_y = mean(sim_y) ; std_y = std(sim_y) ;

% Get target design
obs_x = linspace(0,1,obs_x_size)' * xrange + xmin;

% Emulator mean
mean_sim = @(a,varargin) repmat([0 0],size(a,1),1) ; 

% Define grid of targets
ngrid = 20 ; % Number of points in grid
mincost = 96; maxcost = 352; % Minimum and maximum possible cost
des_obs_grid = [zeros(ngrid,1) linspace(mincost,maxcost,ngrid)'] ;

% Perform CTO for each point in the grid
for ii = 1 : ngrid
    
    % Select the target design to be the point chosen in CTO (elsewhere)
    obs_y = repmat(des_obs_grid(ii,:),obs_x_size,1);
    
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
        'EmulatorCovHypers',EmulatorCovHypers,...
        'des_discrep',false,...
        'obs_var',obs_var,...
        'obs_discrep_use_MLEs',false,...
        'doplot',doplot,...
        'verbose',verbose,...
        'obs_var_est',[true false],...
        'CTO',true);

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
    
    results{ii} = res;
    save temp
    
end

locstr = [dpath,'stored_data\'...
    '2019-11-05_CTO_costgrid'];
% save(locstr,'results');

%% WTA Get posterior predictive results for cost grid estimate of PF
clc ; clearvars -except dpath ; close all ;

% Load the results
locstr = [dpath,'stored_data\'...
    '2019-11-05_CTO_costgrid'];
load(locstr);

% Define inputs mins and ranges 
xmin = results{1}.settings.min_x;
xrange = results{1}.settings.range_x;
tmin = results{1}.settings.min_t1;
trange = results{1}.settings.range_t1;

% Build whatever matrices can be done outside of the main loop
sim_xt = [results{1}.settings.sim_x results{1}.settings.sim_t1 ...
    results{1}.settings.sim_t2];
eta = results{1}.settings.sim_y;
mean_sim = results{1}.settings.mean_sim;
prior_mean_simxt = mean_sim(sim_xt) ;
cov_xt_xt = nan(size(eta,1),size(eta,1),2);
burn_in = results{1}.settings.burn_in;
for ii = 1 : 2
    emulator_rho = results{1}.settings.emulator_rho(:,ii);
    emulator_lambda = results{1}.settings.emulator_lambda(ii);
    cov_xt_xt(:,:,ii) = ...
        gp_cov(emulator_rho, sim_xt, sim_xt, ...
        emulator_lambda,false);
    cov_xt_xt(:,:,ii) = cov_xt_xt(:,:,ii) + ...
        1e-5*eye(size(cov_xt_xt(:,:,ii)));
end

for kk = 1 : length(results)

    % Get posterior samples
    t1 = results{kk}.theta1(:,1); t2 = results{kk}.theta1(:,2);

    % Get set up for loop that will find posterior predictive dist
    pred_xt_all = ([.5*ones(size(t1,1),1) t1 t2] - [0 tmin])./[1 trange];
    prior_mean_xt_all = mean_sim(pred_xt_all) ; 

    % Get posterior predictive distribution for each output
    % Only do maxsimul samples at once, to ease computational expense
    maxsimul = 1000;
    for jj = 1 : ceil(size(t1,1)/maxsimul)
        startpt = (jj-1) * maxsimul + 1 ;
        endpt = min(jj * maxsimul, size(t1,1));
        pred_xt = pred_xt_all(startpt:endpt,:);
        for ii=1:2
            emulator_rho = results{kk}.settings.emulator_rho(:,ii);
            emulator_lambda = results{kk}.settings.emulator_lambda(ii);
            cov_pxt_xt = ...
                gp_cov(emulator_rho, pred_xt, sim_xt, ...
                emulator_lambda,false);
            cov_pxt_pxt = ...
                gp_cov(emulator_rho, pred_xt, pred_xt,...
                emulator_lambda,false);
            pred_mean_xt_current(:,ii) = ...
                prior_mean_xt_all(startpt:endpt,ii) + ...
                cov_pxt_xt * (cov_xt_xt(:,:,ii) \ ...
                (eta(:,ii) - prior_mean_simxt(:,ii))) ; 
            pred_cov_xt_current(:,:,ii) = ...
                cov_pxt_pxt - cov_pxt_xt * ...
                    (cov_xt_xt(:,:,ii) \ cov_pxt_xt') ; 
            fprintf('\n%g of %g, %g of 3 Done\n',jj,...
                ceil(size(t1,1)/maxsimul),ii);
            pred_mean_xt(startpt:endpt,ii) = pred_mean_xt_current(:,ii);
            pred_cov_xt(startpt:endpt,startpt:endpt,ii) = ...
                pred_cov_xt_current(:,:,ii);
        end
    end
    % Put back on original scale
    mean_y = results{kk}.settings.mean_y; 
    std_y = results{kk}.settings.std_y;
    pred_mean_xt_os = pred_mean_xt .* std_y + mean_y;
    mean(pred_mean_xt_os)


    % Now resample from the relevant gaussian dists, for full uncertainty
    for ii = 1:2
        pred_stds_xt(:,ii) = sqrt(diag(pred_cov_xt(:,:,ii)));
        pred_samps_xt(:,ii) = randn(size(pred_mean_xt,1),1) .* ...
            pred_stds_xt(:,ii) + pred_mean_xt(:,ii);
    end
    % Put back on original scale
    pred_samps_xt_os = pred_samps_xt .* std_y + mean_y;
    mean(pred_samps_xt_os)


    %%%%%%%%%%% Make a figure
    % xlims = [0.7135    1.0765 ; 0.6425    1.0275 ; 14.2500   30.7500];
    % ylims = [0   95 ; 0   120 ; 0    3.5000];
    figure('Position',[10 10 600 300],'Color','white');
    for ii = 1:2
        subplot(1,3,ii);
%         histogram(prior_mean_xt_os(:,ii),'Normalization','pdf',...
%             'EdgeColor','none');
%         hold on;
        histogram(pred_mean_xt_os(burn_in:end,ii),'Normalization','pdf',...
            'EdgeColor','none');
        xlabel(sprintf('y_%g',ii));
        set(gca,'Yticklabel',[]);
        set(gca,'Ytick',[]);
    %     xlim(xlims(ii,:));
    %     ylim(ylims(ii,:));
    end

    % Save the posterior predictive output
    model_output.means = pred_mean_xt;
    model_output.vars = cell2mat(...
        arrayfun(@(b)diag(pred_cov_xt(:,:,b)),1:2,'UniformOutput',false));
    results{kk}.model_output = model_output;
    save temp
end

% save(locstr,'results');

%% WTA model validation (via cross validation)
clc ; clearvars -except dpath ; close all ; 

% Set seed for reproducibility
rng(355);

%%% Load raw data
load([dpath,'stored_data\'...
    'raw_dat']);
% Get indices of samples to include in CV
nobs = size(raw_dat,1);
ninc = 20; % number to include
inc_idx = sort(randsample(nobs,ninc));
sim_x = raw_dat(inc_idx,1);
sim_t = raw_dat(inc_idx,2:3);
sim_y = raw_dat(inc_idx,4:6); 
clear raw_dat

%%% Get cross-validation indices
k = 10;
cv_idxs = cvpartition(ninc,'KFold',k);

% Initialize arrays to record results of k-fold CV
RMSE = nan(k,size(sim_y,2));
sim_y_pred = nan(size(sim_y));
sim_y_vars = nan(size(sim_y));

% Perform k-fold CV, get RMSE in each fold
for kk = 1:k
    % Get training and test sets
    sim_x_train = sim_x(cv_idxs.training(kk),:);
    sim_t_train = sim_t(cv_idxs.training(kk),:);
    sim_y_train = sim_y(cv_idxs.training(kk),:);
    sim_x_test = sim_x(cv_idxs.test(kk),:);
    sim_t_test = sim_t(cv_idxs.test(kk),:);
    sim_y_test = sim_y(cv_idxs.test(kk),:);
    
    % Define inputs mins and ranges 
    xmin = min(sim_x_train);
    xrange = range(sim_x_train);
    tmin = min(sim_t_train);
    trange = range(sim_t_train);

    % Get output means and stds
    mean_y = mean(sim_y_train) ; std_y = std(sim_y_train) ;
    
    % Emulator mean
    mean_sim = @(a,varargin) repmat([0 0 0],size(a,1),1) ;	
    
    % Get settings for DCTO, just as an easy way to get GP MLEs
    settings = MCMC_dual_calib_settings(...
        sim_x_train,sim_t_train,[],sim_y_train,...
        [],[],zeros(0,3),...
        [],zeros(0,3),'min_x',xmin,'range_x',xrange,...
        'min_t1',tmin,'range_t1',trange,...
        'min_t2',[],'range_t2',[],...
        'mean_y',mean_y,'std_y',std_y,...
        'obs_discrep',false,...
        'emulator_use',true,...
        'EmulatorMean',mean_sim,...
        'des_discrep',false,...
        'obs_discrep_use_MLEs',false,...
        'CTO',true);
    
    % Get set up for loop that will predict test outputs
    xt_test = ([sim_x_test sim_t_test] - [xmin tmin])./[xrange trange];
    xt_train = ([sim_x_train sim_t_train] - [xmin tmin])./[xrange trange];
    prior_mean_xt_test = settings.mean_sim(xt_test) ; 
    y_train = settings.sim_y;
    prior_mean_xt_train = settings.mean_sim(xt_train) ; 
    cov_train_train = nan(size(xt_train,1),size(xt_train,1),3);
    for jj=1:3
        cov_train_train(:,:,jj) = ...
            gp_cov(settings.emulator_rho(:,jj), xt_train, xt_train,...
            settings.emulator_lambda(jj),false);
        cov_train_train(:,:,jj) = cov_train_train(:,:,jj) + ...
            1e-5*eye(size(cov_train_train(:,:,jj)));
    end
    
    % Get predictive output at each test observation
    for ii=1:3
        cov_test_train = ...
            gp_cov(settings.emulator_rho(:,ii), ...
            xt_test, xt_train, ...
            settings.emulator_lambda(ii),false);
        cov_test_test = ...
            gp_cov(settings.emulator_rho(:,ii), xt_test, xt_test,...
            settings.emulator_lambda(ii),false);
        y_test_pred(:,ii) = ...
            prior_mean_xt_test(:,ii) + ...
            cov_test_train * (cov_train_train(:,:,ii) \ ...
            (y_train(:,ii) - prior_mean_xt_train(:,ii))) ; 
        y_test_cov(:,:,ii) = ...
            cov_test_test - cov_test_train * ...
            (cov_train_train(:,:,ii) \ cov_test_train') ; 
        fprintf('%g of 3 Done\n',ii);
    end
    % Put means back on original scale
    mean_y = settings.mean_y; std_y = settings.std_y;
    y_test_pred_os = y_test_pred .* std_y + mean_y;
    % Get vars, put back on original scale
    y_test_vars = nan( size(y_test_cov,1),size(y_test_cov,3) );
    for ii=1:3
        y_test_vars(:,ii) = diag(y_test_cov(:,:,ii)) ;
    end
    y_test_vars_os = y_test_vars .* std_y.^2;
    
    % Get RMSEs for this fold
    RMSE(kk,:) = sqrt(mean((sim_y_test - y_test_pred_os).^2));
    
    % Record predictions
    sim_y_pred(cv_idxs.test(kk),:) = y_test_pred_os;
    sim_y_vars(cv_idxs.test(kk),:) = y_test_vars_os;
    
    % update
    fprintf(['\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',...
        'STEP %g OF %g DONE\n\n'],kk,k);
end

% Save results
results = struct('RMSEs',RMSE,'predictions',sim_y_pred,...
    'variances',sim_y_vars,'true_vals',sim_y);
locstr = [dpath,'stored_data\'...
    '2020-03-29_GP_validation_CV_results_n',num2str(ninc)];
% save(locstr,'results');

%% WTA get plots for model validation (predicted vs true values)
clc ; clearvars -except dpath ; close all 

% Get names of files storing cross-validation results
files = ls([dpath,'stored_data\2020-03-29_GP_validation_CV_results_n*']);
nfiles = size(files,1);

% Make plot title names
outputs = {'Deflection','Rotation','Cost'};

% Initialize array of figures
figs = gobjects(nfiles,1);

% Loop through all sets of CV results, produce a plot for each
for ii=1:nfiles

    % Initialize figure
    figs(ii) = figure('Position',[40 40 900 250]);

    % Load file
    load(strtrim(files(ii,:)));

    % How big of a data set was used here?
    ninc = size(results.true_vals,1);

    % Get error bar lengths
    err = 1.96 * sqrt(results.variances);

    % Loop through the three outputs
    for kk=1:3
        subplot(1,3,kk);
        errorbar(results.true_vals(:,kk),results.predictions(:,kk),...
            err(:,kk),'*b');
        xlabel(['FE Model: ',outputs{kk}]);
        ylabel(['GP Model: ',outputs{kk}]);
        xlim([min([results.true_vals(:,kk);results.predictions(:,kk)]),...
            max([results.true_vals(:,kk);results.predictions(:,kk)])]);
        ylim([min([results.true_vals(:,kk);results.predictions(:,kk)]),...
            max([results.true_vals(:,kk);results.predictions(:,kk)])]);
        % Add reference diagonal line
        refline(1,0)
    end

    % Add title
    stitle = sprintf('CV results for data set of size %g',ninc);
    sgtitle(stitle);
    
    % White background
    set(figs(ii),'Color','white');
    
    export_fig(stitle);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Everything above this line is actually used to gather results in the   %
% paper. Everything below this line is just workspace used to investigate %
% things.                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get simulation observations from toy ex. as LHC design // not used 
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

%% Gather toy example results on CTO using DCTO code, estimating var
clc ; clearvars -except dpath ; close all ;

% Set model observations
% load([dpath,'Example\Ex_results\'...
% '2019-10-16-raw_dat-3-6-6']);
% load([dpath,'Example\Ex_results\'...
% '2019-10-17-raw_dat-200obs']);
% load([dpath 'Example\Ex_results\'...
% '2019-10-21-raw_dat-30obs'],...
% 'raw_dat');
% load([dpath,'Example\Ex_results\'...
% '2019-10-21-raw_dat-100obs']);
load([dpath,'Example\Ex_results\'...
'2018-05-28-raw_dat-3-12-12']);
% load([dpath,'Example\Ex_results\'...
% '2019-10-24-raw_dat-3-8-8']);

sim_x = raw_dat.sim_xt(:,1);
sim_t = raw_dat.sim_xt(:,2:3);
sim_y = raw_dat.sim_y;
clear raw_dat;

% Settings
obs_var = 0.05; % initial guess of observation error var
obs_x_size  = 3;
doplot = true;
verbose = true;

% Set number of draws, burn_in for each mcmc:
M = 1.2e4; b = .5 ; burn_in = M*b;

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
obs_y = repmat(min(sim_y),obs_x_size,1); % utopia point
% obs_y = repmat([min(sim_y(:,1:2)) 18],obs_x_size,1);

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
    'doplot',doplot,...
    'verbose',verbose,...
    'obs_var_est',true,...
    'CTO',true);

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

% save([dpath,'Example\Ex_results\'...
%     '2019-10-31_CTO_grid3_12_12_varest'],...
%     'res');

%% Look at toy problem results on heatmap including true optimum location
clc ; clearvars -except dpath ; close all;

% load true Pareto front
load([dpath,'Example\Ex_results\'...
    '2018-05-29_true_ctheta-output_nondom'],...
    'ctheta_output_nondom');
% And load results
load([dpath,'Example\Ex_results\'...
    '2019-10-31_CTO_grid3_12_12_varest']);

% Get standardized version
cntrl_mins = res.settings.min_x    ;
cntrl_rngs = res.settings.range_x  ;
calib_mins = res.settings.min_t1   ;
calib_rngs = res.settings.range_t1 ;
output_mns = res.settings.mean_y       ;
output_sds = res.settings.std_y         ;

% Set desired observation
% des_obs_os = [ 0 0 0 ] ; % on original scale
% des_obs_os = obs_y(1,:) ; % on original scale
% des_obs    = (des_obs_os - output_mns)./output_sds;
des_obs = res.settings.obs_y(1,:);
des_obs_os = des_obs .* output_sds + output_mns;

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
set(gcf,'color','white');
xlabel('\theta_1'); ylabel('\theta_2');
% export_fig 'heatmap' -m2
hold on;
burn_in = res.settings.burn_in;
scatter(theta_os(burn_in:end,1),...
    theta_os(burn_in:end,2),5,'.g');
% export_fig 'heatmap_with_results' -m2
% plot(mean(theta_os(burn_in:end,1)),...
%     mean(theta_os(burn_in:end,2)),'.m','MarkerSize',15);   
plot(optim_calib(1),optim_calib(2),'+m','MarkerSize',15,'LineWidth',1.5);
% export_fig 'heatmap_with_results_and_optimum' -m2

%% Gather results on WTA using DCTO code, estimating des_var
clc ; clearvars -except dpath ; close all ; 

%%% Load raw data and desired observation
load([dpath,'stored_data\'...
    'raw_dat']);
sim_x = raw_dat(:,1);
sim_t = raw_dat(:,2:3);
sim_y = raw_dat(:,4:6); 
clear raw_dat

% Settings
obs_var = 0.05; % initial guess of observation error var
obs_x_size  = 3; % Number of target observations
doplot = true;
verbose = true;

% Set number of draws, burn_in for each mcmc:
M = 1e4; b = .5 ; burn_in = M*b;

% Define inputs mins and ranges 
xmin = min(sim_x);
xrange = range(sim_x);
tmin = min(sim_t);
trange = range(sim_t);

% Get output means and stds
mean_y = mean(sim_y) ; std_y = std(sim_y) ;

% Get target design
obs_x = linspace(0,1,obs_x_size)' * xrange + xmin;
% Get target output
load([dpath,'stored_data\'...
    '2018-07-26_elbow_des_obs_d-p2']);
obs_y = repmat(des_obs_new_os,size(obs_x,1),1);
% obs_y = repmat([0 0 0],size(obs_x,1),1); 
% Get est of target error variance
obs_y_std = (obs_y - mean_y)./std_y;
% obs_var = (min(obs_y_std,[],1)-min((sim_y-mean_y)./std_y)).^2;
% obs_y = repmat(min(sim_y),obs_x_size,1); % Utopia point
% obs_y = repmat([min(sim_y(:,1:2)) 18],obs_x_size,1);

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
    'doplot',doplot,...
    'verbose',verbose,...
    'obs_var_est',true,...
    'CTO',true);

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

% save([dpath,'stored_data\'...
%     '2019-10-31_CTO_after_PCTO_desvarest'],...
%     'res');

%% WTA prior predictive distribution vs posterior predictive distribution
clc ; clearvars -except dpath ; close all ;

%%% Load prior predictive results
load([dpath,'stored_data\'...
    '2018-09-03_prior_pred_distrib'],...
    'prior_pred_dist');
prsamps = prior_pred_dist.prior_pred_pts;
clear prior_pred_dist;

%%% Load calib results
% locstr = [dpath,'stored_data\'...
%     '2018-07-27_discrepancy_d-elbow_d-p2'];
locstr = [dpath,'stored_data\'...
    '2019-10-31_CTO_on_utopia_pt_desvarest'];
load(locstr);
posamps = results.model_output.by_sample_est;
des_obs = results.settings.desired_obs;
clear results; 

%%% Make figure using histograms
f=figure('pos',[10 10  640.0000  300]);
% Deflection
subplot(1,3,1);
[p,x]=ksdensity(posamps(:,1));
plot(x,p,'LineWidth',2);
%histogram(posamps(:,1),'Normalization','pdf','Edgecolor','none'); 
hold on;
[p,x]=ksdensity(prsamps(:,1));
plot(x,p,'--','LineWidth',2);
%histogram(prsamps(:,1),'Normalization','pdf','Edgecolor','none');
text(0.05,.9,...
    'Deflection','VerticalAlignment','bottom','Units','normalized');
text(1.715,102,'Rotation','VerticalAlignment','bottom');
xlim([0.6 0.85]);
ylim([0 110]);
line([des_obs(1) des_obs(1)],ylim,'Color','black','Linestyle',':',...
    'linewidth',2);
% Rotation
subplot(1,3,2);
[p,x]=ksdensity(posamps(:,2));
plot(x,p,'LineWidth',2);
%histogram(posamps(:,2),'Normalization','pdf','Edgecolor','none'); 
hold on;
[p,x]=ksdensity(prsamps(:,2));
plot(x,p,'--','LineWidth',2);
%histogram(prsamps(:,2),'Normalization','pdf','Edgecolor','none');
text(0.05,.9,...
    'Rotation','VerticalAlignment','bottom','Units','normalized');
xlim([0.075,0.105])
ylim([0 700]);
line([des_obs(2) des_obs(2)],ylim,'Color','black','Linestyle',':',...
    'linewidth',2);
% Cost
subplot(1,3,3);
[p,x]=ksdensity(posamps(:,3));
plot(x,p,'LineWidth',2);
%histogram(posamps(:,3),'Normalization','pdf','Edgecolor','none'); 
hold on;
[p,x]=ksdensity(prsamps(:,3));
plot(x,p,'--','LineWidth',2);
%histogram(prsamps(:,3),'Normalization','pdf','Edgecolor','none');
text(0.05,.9,...
    'Cost','VerticalAlignment','bottom','Units','normalized');
ylim([0 .0700]);
line([des_obs(3) des_obs(3)],ylim,'Color','black','Linestyle',':',...
    'linewidth',2);
% Add suptitle
% st=suptitle('Prior and posterior predictive distributions');
% st.Position=[0.5 -.1 0];
lg=legend('Posterior','Prior','Target','Location','northeast');
% pos=lg.Position; lg.Position = pos + [0.017 0.065 0 0];
flushLegend(lg,'northeast');

%%% Save
set(f, 'Color','white');
% export_fig FIG_prior_vs_posterior_dist -m3 -painters

%% Get toy example prior vs posterior predictive results
clc ; clearvars -except dpath ; close all ;

% Load results
% locstr = [dpath,'Example\Ex_results\'...
%     '2019-10-24_CTO_with_emulator_grid3_12_12'];
load(locstr);

% Define inputs mins and ranges 
xmin = 1.95;
xrange = .1;
tmin = [0 0];
trange = [3 6];

% Get posterior samples
burn_in = res.settings.burn_in;
t1 = res.theta1(:,1); t2 = res.theta1(:,2);

% Get set up for loop that will find posterior predictive dist
pred_xt = ([2*ones(size(t1,1),1) t1 t2] - [xmin tmin])./[xrange trange];
sim_xt = [res.settings.sim_x res.settings.sim_t1 res.settings.sim_t2];
prior_mean_xt = res.settings.mean_sim(pred_xt) ; 
eta = res.settings.sim_y;
prior_mean_simxt = res.settings.mean_sim(sim_xt) ; 

% Get posterior predictive distribution for each output
for ii=1:3
    cov_pxt_xt = ...
        gp_cov(res.settings.emulator_rho(:,ii), pred_xt, sim_xt, ...
        res.settings.emulator_lambda(ii),false);
    cov_xt_xt = gp_cov(res.settings.emulator_rho(:,ii), sim_xt, sim_xt, ...
        res.settings.emulator_lambda(ii),false);
    cov_xt_xt = cov_xt_xt + 1e-5*eye(size(cov_xt_xt));
    % cov_xt_xt_inv = inv( cov_xt_xt ) ;  
    cov_pxt_pxt = ...
        gp_cov(res.settings.emulator_rho(:,ii), pred_xt, pred_xt,...
        res.settings.emulator_lambda(ii),false);
    pred_mean_xt(:,ii) = ...
        prior_mean_xt(:,ii) + cov_pxt_xt * (cov_xt_xt \ ...
        (eta(:,ii) - prior_mean_simxt(:,ii))) ; 
    pred_cov_xt(:,:,ii) = ...
        cov_pxt_pxt - cov_pxt_xt * (cov_xt_xt \ cov_pxt_xt') ; 
    fprintf('\n%g of 3 Done\n',ii);
end
% Put back on original scale
mean_y = res.settings.mean_y; std_y = res.settings.std_y;
pred_mean_xt_os = pred_mean_xt .* std_y + mean_y;
mean(pred_mean_xt_os)

% Take a look
% histogram(pred_mean_xt_os(:,1),'Normalization','pdf');

% Now resample from the relevant gaussian dists, to get full uncertainty
for ii = 1:3
    pred_stds_xt(:,ii) = sqrt(diag(pred_cov_xt(:,:,ii)));
    pred_samps_xt(:,ii) = randn(size(pred_mean_xt,1),1) .* ...
        pred_stds_xt(:,ii) + pred_mean_xt(:,ii);
end
% Put back on original scale
pred_samps_xt_os = pred_samps_xt .* std_y + mean_y;
mean(pred_samps_xt_os)

% Take a look
% hold on;
% histogram(pred_samps_xt_os(:,1),'Normalization','pdf');


%%%%%%% Get prior predictive output 

% Number of samples
m = 1e4 ; 

% Get theta values
thetavals = rand(m,2) .* trange + tmin ; 
t1 = thetavals(:,1) ; t2 = thetavals(:,2); 

% Get set up for loop that will find prior predictive dist
pred_xt = ([2*ones(size(t1,1),1) t1 t2] - [xmin tmin])./[xrange trange];
prior_mean_xt = res.settings.mean_sim(pred_xt) ; 
clear pred_mean_xt pred_cov_xt pred_stds_xt pred_samps_xt

% Get posterior predictive distribution for each output
for ii=1:3
    cov_pxt_xt = ...
        gp_cov(res.settings.emulator_rho(:,ii), pred_xt, sim_xt, ...
        res.settings.emulator_lambda(ii),false);
    cov_xt_xt = gp_cov(res.settings.emulator_rho(:,ii), sim_xt, sim_xt, ...
        res.settings.emulator_lambda(ii),false);
    cov_xt_xt = cov_xt_xt + 1e-5*eye(size(cov_xt_xt));
    % cov_xt_xt_inv = inv( cov_xt_xt ) ;  
    cov_pxt_pxt = ...
        gp_cov(res.settings.emulator_rho(:,ii), pred_xt, pred_xt,...
        res.settings.emulator_lambda(ii),false);
    pred_mean_xt(:,ii) = ...
        prior_mean_xt(:,ii) + cov_pxt_xt * (cov_xt_xt \ ...
        (eta(:,ii) - prior_mean_simxt(:,ii))) ; 
    pred_cov_xt(:,:,ii) = ...
        cov_pxt_pxt - cov_pxt_xt * (cov_xt_xt \ cov_pxt_xt') ; 
    fprintf('\n%g/3 done\n',ii);
end
% Put back on original scale
mean_y = res.settings.mean_y; std_y = res.settings.std_y;
prior_mean_xt_os = pred_mean_xt .* std_y + mean_y;
mean(prior_mean_xt_os)

% Take a look
% histogram(prior_mean_xt_os(:,1),'Normalization','pdf');

% Now resample from the relevant gaussian dists, to get full uncertainty
for ii = 1:3
    pred_stds_xt(:,ii) = sqrt(diag(pred_cov_xt(:,:,ii)));
    pred_samps_xt(:,ii) = randn(size(pred_mean_xt,1),1) .* ...
        pred_stds_xt(:,ii) + pred_mean_xt(:,ii);
end
% Put back on original scale
prior_samps_xt_os = pred_samps_xt .* std_y + mean_y;
mean(prior_samps_xt_os)

% Take a look
% hold on;
% histogram(prior_samps_xt_os(:,1),'Normalization','pdf');

%%%%%%%%%%% Make a figure
% xlims = [0.7135    1.0765 ; 0.6425    1.0275 ; 14.2500   30.7500];
% ylims = [0   95 ; 0   120 ; 0    3.5000];
fig1 = figure('Position',[10 10 900 300],'Color','white');
for ii = 1:3
    subplot(1,3,ii);
    histogram(prior_mean_xt_os(:,ii),'Normalization','pdf',...
        'EdgeColor','none');
    hold on;
    histogram(pred_mean_xt_os(:,ii),'Normalization','pdf',...
        'EdgeColor','none');
    xlabel(sprintf('y_%g',ii));
    set(gca,'Yticklabel',[]);
    set(gca,'Ytick',[]);
%     xlim(xlims(ii,:));
%     ylim(ylims(ii,:));
end
% export_fig 'posterior_predictive_dists_without_GP_uncertainty' -m2

fig2 = figure('Position',[10 10 900 300],'Color','white');
for ii = 1:3
    subplot(1,3,ii);
    histogram(prior_samps_xt_os(:,ii),'Normalization','pdf',...
        'EdgeColor','none');
    hold on;
    histogram(pred_samps_xt_os(:,ii),'Normalization','pdf',...
        'EdgeColor','none');
    xlabel(sprintf('y_%g',ii));
    set(gca,'Yticklabel',[]);
    set(gca,'Ytick',[]);
%     xlim(xlims(ii,:));
%     ylim(ylims(ii,:));
end
% export_fig 'posterior_predictive_dists_with_GP_uncertainty' -m2

%%%%%%%%%% Make a gif
for jj = 1 : 240
    
    for ii = 1:3
        subplot(1,3,ii);
        histogram(prior_mean_xt_os(:,ii),'Normalization','pdf',...
            'EdgeColor','none');
        hold on;
        histogram(pred_mean_xt_os((jj*25/5):(jj*25),ii),...
            'Normalization','pdf',...
            'EdgeColor','none');
        xlabel(sprintf('y_%g',ii));
        set(gca,'Yticklabel',[]);
        set(gca,'Ytick',[]);
        xlim(xlims(ii,:));
        ylim(ylims(ii,:));
        pause(.01);
        hold off
    end
    
end

%% Get WTA prior vs posterior predictive results
clc ; clearvars -except dpath ; close all ;

% Load results
locstr = [dpath,'stored_data\'...
    '2019-10-31_CTO_on_utopia_pt_desvarest'];
load(locstr);

% Define inputs mins and ranges 
xmin = res.settings.min_x;
xrange = res.settings.range_x;
tmin = res.settings.min_t1;
trange = res.settings.range_t1;
tempval = 310; % We'll look at the outputs at this temperature setting

% Get posterior samples
burn_in = res.settings.burn_in;
t1 = res.theta1(:,1); t2 = res.theta1(:,2);

% Get set up for loop that will find posterior predictive dist
pred_xt = ([tempval * ones(size(t1,1),1) t1 t2] - ...
    [xmin tmin])./[xrange trange];
sim_xt = [res.settings.sim_x res.settings.sim_t1 res.settings.sim_t2];
prior_mean_xt = res.settings.mean_sim(pred_xt) ; 
eta = res.settings.sim_y;
prior_mean_simxt = res.settings.mean_sim(sim_xt) ; 

% Get posterior predictive distribution for each output
for ii=1:3
    cov_pxt_xt = ...
        gp_cov(res.settings.emulator_rho(:,ii), pred_xt, sim_xt, ...
        res.settings.emulator_lambda(ii),false);
    cov_xt_xt = gp_cov(res.settings.emulator_rho(:,ii), sim_xt, sim_xt, ...
        res.settings.emulator_lambda(ii),false);
    cov_xt_xt = cov_xt_xt + 1e-5*eye(size(cov_xt_xt));
    % cov_xt_xt_inv = inv( cov_xt_xt ) ;  
    cov_pxt_pxt = ...
        gp_cov(res.settings.emulator_rho(:,ii), pred_xt, pred_xt,...
        res.settings.emulator_lambda(ii),false);
    pred_mean_xt(:,ii) = ...
        prior_mean_xt(:,ii) + cov_pxt_xt * (cov_xt_xt \ ...
        (eta(:,ii) - prior_mean_simxt(:,ii))) ; 
    pred_cov_xt(:,:,ii) = ...
        cov_pxt_pxt - cov_pxt_xt * (cov_xt_xt \ cov_pxt_xt') ; 
    fprintf('\nPosterior predictions %g of 3 Done\n',ii);
end
% Put back on original scale
mean_y = res.settings.mean_y; std_y = res.settings.std_y;
pred_mean_xt_os = pred_mean_xt .* std_y + mean_y;
mean(pred_mean_xt_os)

% Take a look
histogram(pred_mean_xt_os(burn_in:end,3),'Normalization','pdf');

% Now resample from the relevant gaussian dists, to get full uncertainty
for ii = 1:3
    pred_stds_xt(:,ii) = sqrt(diag(pred_cov_xt(:,:,ii)));
    pred_samps_xt(:,ii) = randn(size(pred_mean_xt,1),1) .* ...
        pred_stds_xt(:,ii) + pred_mean_xt(:,ii);
end
% Put back on original scale
pred_samps_xt_os = pred_samps_xt .* std_y + mean_y;
mean(pred_samps_xt_os)

% Take a look
hold on;
histogram(pred_samps_xt_os(burn_in:end,3),'Normalization','pdf');


%%%%%%% Get prior predictive output 

% Number of samples
m = 1e4 ; 

% Get theta values
thetavals = rand(m,2) .* trange + tmin ; 
t1 = thetavals(:,1) ; t2 = thetavals(:,2); 

% Get set up for loop that will find prior predictive dist
pred_xt = ([tempval * ones(size(t1,1),1) t1 t2] - ...
    [xmin tmin])./[xrange trange];
prior_mean_xt = res.settings.mean_sim(pred_xt) ; 
clear pred_mean_xt pred_cov_xt pred_stds_xt pred_samps_xt

% Get posterior predictive distribution for each output
for ii=1:3
    cov_pxt_xt = ...
        gp_cov(res.settings.emulator_rho(:,ii), pred_xt, sim_xt, ...
        res.settings.emulator_lambda(ii),false);
    cov_xt_xt = gp_cov(res.settings.emulator_rho(:,ii), sim_xt, sim_xt, ...
        res.settings.emulator_lambda(ii),false);
    cov_xt_xt = cov_xt_xt + 1e-5*eye(size(cov_xt_xt));
    % cov_xt_xt_inv = inv( cov_xt_xt ) ;  
    cov_pxt_pxt = ...
        gp_cov(res.settings.emulator_rho(:,ii), pred_xt, pred_xt,...
        res.settings.emulator_lambda(ii),false);
    pred_mean_xt(:,ii) = ...
        prior_mean_xt(:,ii) + cov_pxt_xt * (cov_xt_xt \ ...
        (eta(:,ii) - prior_mean_simxt(:,ii))) ; 
    pred_cov_xt(:,:,ii) = ...
        cov_pxt_pxt - cov_pxt_xt * (cov_xt_xt \ cov_pxt_xt') ; 
    fprintf('\nPrior predictions %g of 3 done\n',ii);
end
% Put back on original scale
mean_y = res.settings.mean_y; std_y = res.settings.std_y;
prior_mean_xt_os = pred_mean_xt .* std_y + mean_y;
mean(prior_mean_xt_os)

% Take a look
% histogram(prior_mean_xt_os(:,1),'Normalization','pdf');

% Now resample from the relevant gaussian dists, to get full uncertainty
for ii = 1:3
    pred_stds_xt(:,ii) = sqrt(diag(pred_cov_xt(:,:,ii)));
    pred_samps_xt(:,ii) = randn(size(pred_mean_xt,1),1) .* ...
        pred_stds_xt(:,ii) + pred_mean_xt(:,ii);
end
% Put back on original scale
prior_samps_xt_os = pred_samps_xt .* std_y + mean_y;
mean(prior_samps_xt_os)

% Take a look
% hold on;
% histogram(prior_samps_xt_os(:,1),'Normalization','pdf');

%%%%%%%%%%% Make a figure
% xlims = [0.7135    1.0765 ; 0.6425    1.0275 ; 14.2500   30.7500];
% ylims = [0   95 ; 0   120 ; 0    3.5000];
fig1 = figure('Position',[10 10 900 300],'Color','white');
for ii = 1:3
    subplot(1,3,ii);
    histogram(prior_mean_xt_os(:,ii),'Normalization','pdf',...
        'EdgeColor','none');
    hold on;
    histogram(pred_mean_xt_os(:,ii),'Normalization','pdf',...
        'EdgeColor','none');
    xlabel(sprintf('y_%g',ii));
    set(gca,'Yticklabel',[]);
    set(gca,'Ytick',[]);
%     xlim(xlims(ii,:));
%     ylim(ylims(ii,:));
end
% export_fig 'posterior_predictive_dists_without_GP_uncertainty' -m2

fig2 = figure('Position',[10 10 900 300],'Color','white');
for ii = 1:3
    subplot(1,3,ii);
    histogram(prior_samps_xt_os(:,ii),'Normalization','pdf',...
        'EdgeColor','none');
    hold on;
    histogram(pred_samps_xt_os(:,ii),'Normalization','pdf',...
        'EdgeColor','none');
    xlabel(sprintf('y_%g',ii));
    set(gca,'Yticklabel',[]);
    set(gca,'Ytick',[]);
%     xlim(xlims(ii,:));
%     ylim(ylims(ii,:));
end
% export_fig 'posterior_predictive_dists_with_GP_uncertainty' -m2

%%%%%%%%%% Make a gif
for jj = 1 : 240
    
    for ii = 1:3
        subplot(1,3,ii);
        histogram(prior_mean_xt_os(:,ii),'Normalization','pdf',...
            'EdgeColor','none');
        hold on;
        histogram(pred_mean_xt_os((jj*25/5):(jj*25),ii),...
            'Normalization','pdf',...
            'EdgeColor','none');
        xlabel(sprintf('y_%g',ii));
        set(gca,'Yticklabel',[]);
        set(gca,'Ytick',[]);
%         xlim(xlims(ii,:));
%         ylim(ylims(ii,:));
        pause(.01);
        hold off
    end
    
end
