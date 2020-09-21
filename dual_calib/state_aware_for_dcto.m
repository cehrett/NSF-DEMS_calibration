% State-aware for DCTO
% This file is for gathering state-aware results to use in the DCTO paper.

%% Set path string and add paths
clc; clear all; close all;

direc = pwd; 
if direc(1)=='C' 
    dpath = 'C:\Users\Carl\Documents\MATLAB\NSF_DEMS\NSF-DEMS_calibration\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc;

% Add paths
addpath(dpath);
addpath([dpath,'stored_data']);
addpath([dpath,'Example']);
addpath([dpath,'state-aware']);
addpath([dpath,'dual_calib']);

% Change dir
cd(dpath);

% Confirm
disp(dpath)

%% Perform state-aware calibration on system from DCTO paper without x
clc ; clearvars -except dpath ; close all ;
% This is a messy combination of two very different chunks of code.
% Probably contains lots of things that ought to be deleted. Not sure if
% the logic here does what it appears to do. Tread with caution.

% Set discrepancy version
discrep=0;

% Set initial and final observation set sizes
obs_initial_size = 10;
obs_final_size = 20;

% Settings
modular = true;
informative_targets = false;
des_discrep = false;
obs_discrep_use_MLEs = false;
obs_var = 0.05; % observation error
des_x_size = 15;
doplot = true;
verbose = false;
obs_var_est = false;

% Set number of draws, burn_in for each mcmc:
M = 8e3; b = .5 ; burn_in = M*b;

% Define inputs mins and ranges 
xmin = .5;
xrange = .5;
t1min = 1.5;
t1range = 3;
t2min = 0;
t2range = 5;

% Observation and target discrepancy covariance hyperparameters
obs_rho_beta_params = [8,1];
obs_lambda_gam_params = [8,4];
des_rho_beta_params = [2,.4];
des_lambda_gam_params = [40,40];

% True theta1 function
high_theta1 = 2.25 ; low_theta1 = 1.5 ; 
%         t1_fn = @(x) high_theta1 - ...
%             (high_theta1-low_theta1) * ...
%             exp(40*((x-t2min)/t2range)-20)./...
%             (1+exp(40*((x-t2min)/t2range)-20));
t1_fn = @(x) high_theta1 - ...
    (high_theta1-low_theta1) * ...
    exp(40*((x-t2min)/t2range)-20)./...
    (1+exp(40*((x-t2min)/t2range)-20));


% Define comp. model & truth (version w/ nonstandardized output)
model_fn_ns = @(x_t1_t2) dual_calib_example_fn(...
    x_t1_t2(:,1),xmin,xrange,x_t1_t2(:,2),t1min,t1range,x_t1_t2(:,3),t2min,t2range,...
    0,1,0,true) ; 
true_phenomenon_ns = @(x,t2) dual_calib_example_fn(...
    x,xmin,xrange,...
    (t1_fn(t2*t2range+t2min)-t1min)./t1range,t1min,t1range,...
    t2,t2min,t2range,0,1,discrep,true);


% Get mean and std using comp. model, define standardized version
X=lhsdesign(1000,3); 
Y=model_fn_ns([X(:,1) X(:,2) X(:,3)]);
mean_y = mean(Y) ; std_y = std(Y) ;
model_fn = @(x_t1_t2) (model_fn_ns(x_t1_t2) - mean_y)./std_y;
true_phenomenon = ...
    @(x,t2) (true_phenomenon_ns(x,t2) - mean_y)./std_y;

% Get initial design
X = lhsdesign(obs_initial_size,2);
obs_x = X(:,1) * xrange + xmin ; 
obs_t2 = X(:,2) * t2range + t2min ;

% Get "real" observations without noise but with discrepancy
obs_y_noiseless = true_phenomenon_ns((obs_x-xmin)./xrange,...
    (obs_t2-t2min)./t2range);

% Now add N(0,obs_var) noise to STANDARDIZED observations
obs_y = obs_y_noiseless+ ...
    randn(obs_initial_size,1) * sqrt(obs_var) * std_y;

% Now set desired observations
des_x = linspace(0,1,des_x_size)' * xrange + xmin;
if informative_targets    
    % Find true minimum at des_x
    fn=@(t2)dual_calib_example_fn(.5,xmin,xrange,...
        t1_fn(t2),t1min,t1range,...
        t2,t2min,t2range,0,1,discrep,false);
    theta2 = ...
        fmincon(fn,rand*t2range + t2min,...
        [],[],[],[],t2min,t2min+t2range);
        true_min = dual_calib_example_fn(des_x,xmin,xrange,...
            t1_fn(theta2) * ones(size(des_x)),t1min,t1range,...
            theta2,t2min,t2range,0,1,discrep,false);
        des_y = true_min - 2*sqrt(obs_var) ; 
else
    des_y = 0 * ones(size(des_x,1),1);
    % Get est of target error variance
    des_y_std = (des_y - mean_y)/std_y;
    des_var = (min(des_y_std)-min((Y-mean_y)/std_y))^2;
end

% Since we are not using emulator, empty out simulator observations
sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

% And get the average discrepancy value to use as the mean
int_fn =@(x,t) true_phenomenon(x,t) - ...
    model_fn([x (t1_fn(t*t2range+t2min)-t1min)./t1range t]) ;
avg_disc = integral2(int_fn,0,1,0,1) ; 
fprintf('Average observation discrepancy: %g\n',avg_disc);
obs_discrep_mean = @(x,t) avg_disc * ones(size(x)) ; 

% Emulator mean
mean_sim = @(abc) model_fn(abc) ; 

% Set inputs for calibration
desired_obs = obs_y ; 
output_mean = mean_y;
output_sd   = std_y;
lambda_delta = 1/(.5)^2;
cntrl_input = obs_t2;

eta = @(xx,theta1,~,...
    input_cntrl_min,input_cntrl_range,...
    input_calib_min,input_calib_range,...
    output_mean,output_sd,~) dual_calib_example_fn(...
    0.75,.5,.5,...
    theta1,input_calib_min,input_calib_range,...
    xx,input_cntrl_min,input_cntrl_range,...
    output_mean,output_sd,discrep,true);

% Get settings
settings = MCMC_sa_settings(desired_obs,lambda_delta,cntrl_input,...
    'input_cntrl_min',t2min,...
    'input_cntrl_range',t2range,...
    'which_sa',1,...
    'dim_theta1',1,...
    'dim_theta2',0,...
    'input_calib_min',t1min,...
    'input_calib_range',t1range,...
    'eta',eta,...
    'output_mean',output_mean,...
    'output_sd',output_sd,...
    'nu_theta_prior_param',100,...
    'nu_delta_prior_param',100,...
    'lambda_theta_hypers',[1, 1],...
    'M',2e4,...
    'burn_in',1/4);

% Perform calibration
results = MCMC_state_aware_discrep_true_fn ( settings ) ;

% Compare results to truth
res=mean(results.theta1(10001:end,:)) * results.settings.input_calib_range;
truth=t1_fn(cntrl_input);
figure();
plot(cntrl_input,res,'bo'); hold on ; plot(cntrl_input,truth,'ro') ; 

% Save results
% save([dpath,'state-aware\state-aware_results\'...
%     '2018-11-29_sa_univariate_system'],...
%     'results');


%% Examine marginal posteriors of theta1 at grid of control settings
burn_in = results.settings.burn_in;
cal_rng = results.settings.input_calib_range;
cal_min = results.settings.input_calib_min;

samps = results.theta1(burn_in:end,:) * cal_rng + cal_min;
meantheta = mean(samps);
truetheta = t1_fn(results.settings.cntrl_input);
for ii=1:size(samps,2)
    subplot(2,4,ii);
    ksdensity(samps(:,ii));
    xlim([0.5,4]);ylim([0,1.75]);
    hold on;
    plot([meantheta(ii) meantheta(ii)],get(gca,'YLim'));
    plot([truetheta(ii) truetheta(ii)],get(gca,'YLim'));
end

%% Perform state-aware calibration on system from DCTO paper w/ x included
clc ; clearvars -except dpath ; close all ;
% This is a messy combination of two very different chunks of code.
% Probably contains lots of things that ought to be deleted. Not sure if
% the logic here does what it appears to do. Tread with caution.

% Set discrepancy version
discrep=0;

% Set initial and final observation set sizes
obs_initial_size = 10;
obs_final_size = 20;

% Settings
obs_var = 0.05; % observation error
des_x_size = 15;
doplot = true;
verbose = false;

% Set number of draws, burn_in for each mcmc:
M = 3e4; b = .5 ; burn_in = M*b;

% Define inputs mins and ranges 
xmin = .5;
xrange = .5;
t1min = 1.5
t1range = 3
t2min = 0;
t2range = 5;

% True theta1 function
high_theta1 = 2.5 ; low_theta1 = 1.75 ; 
t1_fn = @(x) high_theta1 - ...
    (high_theta1-low_theta1) * ...
    exp(40*((x-t2min)/t2range)-20)./...
    (1+exp(40*((x-t2min)/t2range)-20));

% Define comp. model & truth (version w/ nonstandardized output)
% Note these take inputs on normalized scale
model_fn_ns = @(x_t1_t2) dual_calib_example_fn(...
    x_t1_t2(:,1:(size(x_t1_t2,2)/3)),xmin,xrange,...
    x_t1_t2(:,(size(x_t1_t2,2)/3)+1:(2*size(x_t1_t2,2)/3)),t1min,t1range,...
    x_t1_t2(:,(2*size(x_t1_t2,2)/3)+1:(size(x_t1_t2,2))),t2min,t2range,...
    0,1,... %output mean and std to use for rescaling (ie don't rescale)
    0,... % discrep
    true ) ; % rescale inputs using min and range?
true_phenomenon_ns = @(x,t2) dual_calib_example_fn(... %takes inputs on stdzd scale, returns on original scale 
    x,xmin,xrange,...
    (t1_fn(t2*t2range+t2min)-t1min)./t1range,t1min,t1range,...
    t2,t2min,t2range,...
    0,1,... %output mean and std to use for rescaling (ie don't rescale)
    discrep,...
    true);

% Get mean and std using comp. model, define standardized version
% Note these take inputs on normalized scale
X=lhsdesign(1000,3); 
Y=model_fn_ns([X(:,1) X(:,2) X(:,3)]);
mean_y = mean(Y) ; std_y = std(Y) ;
model_fn = @(x_t1_t2) (model_fn_ns(x_t1_t2) - mean_y)./std_y;
true_phenomenon = ...
    @(x,t2) (true_phenomenon_ns(x,t2) - mean_y)./std_y;

% Get initial design
X = lhsdesign(obs_initial_size,2);
obs_x = X(:,1) * xrange + xmin ; 
obs_t2 = X(:,2) * t2range + t2min ;

% Get "real" observations without noise but with discrepancy
obs_y_noiseless = true_phenomenon_ns((obs_x-xmin)./xrange,...
    (obs_t2-t2min)./t2range); % Note inputs are convered to normalizd scale

%%%  TEMPORARILY REMOVED %%%
% Now add N(0,obs_var) noise to STANDARDIZED observations. Here that noise
% is rescaled and added to the non-standardized observations.
% obs_y = obs_y_noiseless+ ...
%     randn(obs_initial_size,1) * sqrt(obs_var) * std_y;
obs_y = obs_y_noiseless

% Since we are not using emulator, empty out simulator observations
sim_x = [] ; sim_t1 = [] ; sim_t2 = [] ; sim_y = [] ;

% And get the average discrepancy value to use as the mean
int_fn =@(x,t) true_phenomenon(x,t) - ...
    model_fn([x (t1_fn(t*t2range+t2min)-t1min)./t1range t]) ;
avg_disc = integral2(int_fn,0,1,0,1) ; 
fprintf('Average observation discrepancy: %g\n',avg_disc);
obs_discrep_mean = @(x,t) avg_disc * ones(size(x)) ; 

% Emulator mean
mean_sim = @(abc) model_fn(abc) ; 

% Set inputs for calibration
output_mean = mean_y;
output_sd   = std_y;
% lambda_delta = 1/(.5)^2; % Don't know whether this val is appropriate
lambda_delta = 1/(.05)^2; % Don't know whether this val is appropriate
cntrl_input = [ obs_x obs_t2 ];

eta = @(xx,theta1,~,...
    input_cntrl_min,input_cntrl_range,...
    input_calib_min,input_calib_range,...
    output_mean,output_sd,~) dual_calib_example_fn(...
    xx(:,1),input_cntrl_min(1),input_cntrl_range(1),...
    theta1,input_calib_min,input_calib_range,...
    xx(:,2),input_cntrl_min(2),input_cntrl_range(2),...
    output_mean,output_sd,discrep,true);

% Get settings
settings = MCMC_sa_settings(obs_y,lambda_delta,cntrl_input,...
    'input_cntrl_min',[xmin t2min],...
    'input_cntrl_range',[xrange t2range],...
    'which_sa',1,...
    'dim_theta1',1,...
    'dim_theta2',0,...
    'input_calib_min',t1min,...
    'input_calib_range',t1range,...
    'eta',eta,...
    'output_mean',output_mean,...
    'output_sd',output_sd,...
    'nu_theta_prior_param',5,...
    'nu_delta_prior_param',10,...
    'lambda_theta_hypers',[50, 50],...
    'M',M,...
    'burn_in',b,...
    'mu_theta',.1);

% Perform calibration
results = MCMC_state_aware_discrep_true_fn ( settings ) ;

% Compare results to truth
res=mean(results.theta1(burn_in:end,:)) * results.settings.input_calib_range + results.settings.input_calib_min;
truth=t1_fn(cntrl_input(:,2));
figure();
plot(cntrl_input(:,2),res,'bo'); hold on ; plot(cntrl_input(:,2),truth,'ro') ; 

% Save results
% save([dpath,'state-aware\state-aware_results\'...
%     '2018-11-29_sa_univariate_system'],...
%     'results');