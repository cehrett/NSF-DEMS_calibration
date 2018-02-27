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
p = 2 ; % Number of chains at each desired obs. to run in parallel
results_par = cell(p,1);
parpool(p);
n=2;
M=20;
all_results = cell(n,1);
design = lhsdesign(n,2) .* [ .74 224 ] ;
% Note: the range of the design was chosen by selecting 0 as the minimum
% values and for the maximum values taking the midpoint of the plausible
% ranges supplied by Evan for deflection and cost.

for jj=1:n

    parfor ii=1:p
        desired_obs = design(jj,:);
        settings = MCMC_settings(M,2,desired_obs);
        [samples,sigma2_rec,Sigma] = MCMC_sigma_prior_joint_prop(settings);
        post_mean_out = em_out(samples,settings.burn_in,settings.obs_x,...
            settings.sim_xt,settings.eta,settings.output_sds,...
            settings.output_means,settings.omega,settings.rho,...
            settings.lambda);
        results = struct('samples',samples,'sigma2',sigma2_rec,...
            'Sigma',Sigma,'init',settings.init_theta,...
            'desired_data',desired_obs,...
            'sigma2_prior',settings.log_sigma2_prior,...
            'omega_rho_lambda',[settings.omega settings.rho settings.lambda],...
            'proposal',settings.proposal,'nugsize',settings.nugsize,...
            'post_mean_theta',mean(samples(settings.burn_in:end,:)),...
            'post_mean_sigma2',mean(sigma2_rec(settings.burn_in:end,:)),...
            'post_mean_out',post_mean_out);
        results_par{ii} = results;
    end
    
    all_results{jj} = results_par
    
end

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
