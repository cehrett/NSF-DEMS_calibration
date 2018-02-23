% Master file for looped runs

clc; clear all; close all;

%% Add paths
addpath('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration');
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data');
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\');


%% User defined values
M = 1e4;
desired_obs = [.65 0.077 96];
desired_obs = [.65 96];

%% Settings
settings = MCMC_settings (M,3,desired_obs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Joint proposal for theta, prior put on obs variance
[samples,sigma2_rec,Sigma] = MCMC_sigma_prior(M,burn_in,sim_xt,eta,...
    obs_x,y,sigma2,log_sigma2_prior,out_of_range,init_theta,omega,...
    rho,lambda,proposal,nugsize);

results = struct('samples',samples,'sigma2',sigma2_rec,'Sigma',Sigma,...
'init',init_theta,'desired_data',desired_obs,...
'sigma2_prior',log_sigma2_prior,...
'omega_rho_lambda',[omega rho lambda],'proposal',proposal,...
'nugsize',nugsize);

save(['.\NSF DEMS\Phase 1\'...
'results_z0_univObsPrior1'],...
'results');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Joint prop for theta, joint prop for obs var, prior on obs var
[samples,sigma2_rec,Sigma] = MCMC_sigma_prior_joint_prop(settings);

results = struct('samples',samples,'sigma2',sigma2_rec,'Sigma',Sigma,...
'init',init_theta,'desired_data',desired_obs,...
'sigma2_prior',log_sigma2_prior,...
'omega_rho_lambda',[omega rho lambda],'proposal',proposal,...
'nugsize',nugsize);

save(['E:\Carl\Documents\MATLAB\NSF-DEMS_data\'...
'results_z1_hetskedObsPrior1'],...
'results');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parallel processing loop
%% Joint prop for theta, obs var; prior on obs var
% Cell to store data
results_par = cell(3,1);
parpool(3);
parfor ii=1:3
    M = 1e2;
    desired_obs = [.65 0.077 96];
    settings = MCMC_settings (M,desired_obs);
    [samples,sigma2_rec,Sigma] = MCMC_sigma_prior_joint_prop(settings);
    results = struct('samples',samples,'sigma2',sigma2_rec,'Sigma',Sigma,...
        'init',settings.init_theta,'desired_data',desired_obs,...
        'sigma2_prior',settings.log_sigma2_prior,...
        'omega_rho_lambda',[settings.omega settings.rho settings.lambda],...
        'proposal',settings.proposal,'nugsize',settings.nugsize);
    results_par{ii} = results;
end
