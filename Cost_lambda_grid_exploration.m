% Workspace for grid investigation of Cost_lambda

clc; clear all; close all;

%% Add paths
addpath('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration');
addpath('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\stored_data');
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data');
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\');

%% User defined values
m = 2; % Grid size
M = 5e2; % MCMC length
desired_obs = [ 0 0 ];
which_outputs = [ 1 1 0 ] ; %Which of defl, rot, cost

%% Settings
settings = MCMC_settings (M,desired_obs,which_outputs);
settings.doplot = false;
settings.doplot = true;
settings.Cost_lambda = 0; % remove prior on vf, thk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loop through grid of Cost_lambda vals
%gridvals = exp(linspace(0,log(100),m));
gridvals = linspace(0,100,m); % Get grid values
results = cell(m,1); % This will store results
for ii = 1:m
    
    settings.Cost_lambda = gridvals(ii);
    
    % Joint prop for theta, joint prop for obs var, prior on obs var
    [samples,sigma2_rec,Sigma] = MCMC_sigma_prior_joint_prop(settings);
    
    % Collect results
    % trivar_output_settings is useful for getting output at post mean
    trivar_output_settings = MCMC_settings(M,[0 0 0],[1 1 1]);
    post_mean_out = em_out(samples,settings);
    result = struct('samples',samples,'sigma2',sigma2_rec,'Sigma',Sigma,...
        'init',samples(1,:),'desired_obs',desired_obs,...
        'sigma2_prior',settings.log_sigma2_prior,...
        'omega_rho_lambda',[settings.omega settings.rho settings.lambda],...
        'proposal',settings.proposal,'nugsize',settings.nugsize,...
        'post_mean_theta',mean(samples(settings.burn_in:end,:)),...
        'post_mean_sigma2',mean(sigma2_rec(settings.burn_in:end,:)),...
        'post_mean_out',post_mean_out,'which_outputs',which_outputs,...
        'Cost_lambda',settings.Cost_lambda,...
        'log_theta_prior',settings.log_theta_prior,...
        'log_sig_mh_correction',settings.log_sig_mh_correction,...
        'log_mh_correction',settings.log_mh_correction,...
        'burn_in',settings.burn_in);
    
    results{ii} = result;
end

% save(['C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data\'...
% 'results_Cost_lambda_grid_exploration'],...
% 'results');

save(['E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\stored_data\'...
'results_Cost_lambda_grid_exploration'],...
'results');

