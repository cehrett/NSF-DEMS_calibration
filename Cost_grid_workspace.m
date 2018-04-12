% Cost grid pareto band workspace
% This is a workspace in which to make a grid of target costs which will be
% specified with extremely low observation variance, to force the MCMC to
% stay at or very near those cost targets. The performance of the samples
% at each grid point will then be plotted along with error bars for 2 sd.

clc; clear all; close all;

% Set path string
direc = pwd; if direc(1)=='C' 
    dpath = 'C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\';
else
    dpath = 'E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\';
end
clear direc;

%% Add paths
addpath(dpath);
addpath([dpath,'stored_data']);

%% User defined values
m = 10; % Grid size
M = 6e3; % MCMC length
desired_def_rot = [ 0 0 ]; % desired deflection and rotation
which_outputs = [ 1 1 1 ] ; % Which of defl, rot, cost
cost_var = 0.05; % Sets constant obs var for cost

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loop through grid of Costs
gridvals = linspace(96,350,m); % Get grid values
results = cell(m,1); % This will store results
for ii = 1:m
    
    % Set desired_obs for this loop
    desired_obs = [desired_def_rot gridvals(ii)];
    
    % Get settings for MCMC
    settings = MCMC_settings(M,desired_obs,which_outputs);
    
    % Modify settings
    % Many settings must be changed to accommodate the fact that we are
    % treating the observation variance for cost as known.
    % First, modify the proposal density for obs var.
    settings.proposal.sigma2_prop_density = @(x,s) ...
        [exp(mvnrnd(log([x(1) x(2)]),s(1:2,1:2))) x(3)];
    % Now modify the prior for obs var.
    settings.log_sigma2_prior = @(sigma2) -log(prod(sigma2(1:2)));
    % Now modify the initial obs var.
    settings.sigma2 = [settings.sigma2(1:2) cost_var];
    % Now modify the Metrop.-Hastings correction for drawing obs var.
    settings.log_sig_mh_correction = @(sig_s,sig) ...
        log(prod(sig_s(1:2)))-log(prod(sig(1:2)));
    % We don't want an informative prior on VF, thk., so remove that:
    settings.log_theta_prior = @(theta,Cost_lambda) 0 ;
    % Okay, now tell it we want progress plots during the MCMC
    settings.doplot = true;
    % And set the burn_in to what we want it to be
    settings.burn_in=2000;
    
    % Run the MCMC
    [samples,sigma2_rec,Sigma] = MCMC_sigma_prior_joint_prop(settings);
    
    % Collect results
    post_mean_out = em_out(samples,settings);
    result = struct('samples',samples,...
        'sigma2',sigma2_rec,...
        'Sigma',Sigma,...
        'desired_obs',desired_obs,...
        'post_mean_theta',mean(samples(settings.burn_in:end,:)),...
        'post_mean_sigma2',mean(sigma2_rec(settings.burn_in:end,:)),...
        'post_mean_out',post_mean_out,...
        'settings',settings);
    
    results{ii} = result;
    
    save([dpath,'stored_data\'...
        'results_Cost_grid_exploration'],...
        'results');
    
end

save([dpath,'stored_data\'...
    'results_Cost_grid_exploration'],...
    'results');