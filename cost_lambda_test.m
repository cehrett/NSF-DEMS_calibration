clc; clear all; close all;

%% Add paths
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data');
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\');
addpath('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration');
addpath('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\stored_data');


%% User defined values
M = 1e4;
desired_obs = [.65 0.77];
which_outputs = [ 1 1 0 ] ; %Which of defl, rot, cost

%% Settings
settings = MCMC_settings (M,desired_obs,which_outputs);
settings.doplot = false;
settings.doplot = true;
results = cell(3,1);

for ii = [1 10 100]
    settings.Cost_lambda=ii;
    
    [samples,sigma2_rec,Sigma] = MCMC_sigma_prior_joint_prop(settings);


    post_mean_out = em_out(samples,settings.burn_in,settings.obs_x,...
        settings.sim_xt,settings.eta,settings.output_sds,...
        settings.output_means,settings.omega,settings.rho,...
        settings.lambda,which_outputs);
    result = struct('samples',samples,'sigma2',sigma2_rec,'Sigma',Sigma,...
        'init',samples(1,:),'desired_data',desired_obs,...
        'sigma2_prior',settings.log_sigma2_prior,...
        'omega_rho_lambda',[settings.omega settings.rho settings.lambda],...
        'proposal',settings.proposal,'nugsize',settings.nugsize,...
        'post_mean_theta',mean(samples(settings.burn_in:end,:)),...
        'post_mean_sigma2',mean(sigma2_rec(settings.burn_in:end,:)),...
        'post_mean_out',post_mean_out,...
        'Cost_lambda',settings.Cost_lambda);
    
    results{ii}=result;
end

save(['E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\stored_data\'...
'results_Cost_lambda_test'],...
'results');




