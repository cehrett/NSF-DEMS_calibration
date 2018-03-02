% Master file for looped runs

clc; clear all; close all;

%% Add paths
addpath('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration');
addpath('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\stored_data');
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\stored_data');
addpath('C:\Users\carle\Documents\MATLAB\NSF DEMS\Phase 1\');


%% User defined values
M = 1e4;
desired_obs = [.65 0.077 96];
desired_obs = [.65 96];

%% Settings
settings = MCMC_settings (M,desired_obs);

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
n=10;
design = lhsdesign(n,2) .* [ .74 224 ] ;
p = 3 ; % Number of chains at each desired obs. to run in parallel
% Cell to store data
results_par = cell(p,1);
parpool(p);
M=1e4;
all_results = cell(n,1);
% Begin from here if interrupted
load('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\stored_data\30_MCMCs');

% Note: the range of the design was chosen by selecting 0 as the minimum
% values and for the maximum values taking the midpoint of the plausible
% ranges supplied by Evan for deflection and cost.

for jj=4:4

    parfor ii=1:p
        desired_obs = design(jj,:);
        settings = MCMC_settings(M,desired_obs);
        [samples,sigma2_rec,Sigma] = MCMC_sigma_prior_joint_prop(settings);
        post_mean_out = em_out(samples,settings.burn_in,settings.obs_x,...
            settings.sim_xt,settings.eta,settings.output_sds,...
            settings.output_means,settings.omega,settings.rho,...
            settings.lambda);
        results = struct('samples',samples,'sigma2',sigma2_rec,...
            'Sigma',Sigma,'init',settings.init_theta,...
            'desired_data',desired_obs,...
            'sigma2_prior',settings.log_sigma2_prior,...
            'omega_rho_lambda',...
            [settings.omega settings.rho settings.lambda],...
            'proposal',settings.proposal,'nugsize',settings.nugsize,...
            'post_mean_theta',mean(samples(settings.burn_in:end,:)),...
            'post_mean_sigma2',mean(sigma2_rec(settings.burn_in:end,:)),...
            'post_mean_out',post_mean_out);
        results_par{ii} = results;
    end
    
    all_results{jj} = results_par;
    
    fprintf('COMPLETED LOOP %g/%g',jj,n)
    
end

save('E:\Carl\Documents\MATLAB\NSF-DEMS_calibration\stored_data\30_MCMCs');

% Take a look at results
for jj=1:n
    for ii=1:p
        subplot(2,2,1);
        plot(all_results{jj}{ii}.samples(burn_in:end,1),'ko');
        subplot(2,2,2);
        plot(all_results{jj}{ii}.samples(burn_in:end,2),'ko');
        subplot(2,2,3);
        plot(all_results{jj}{ii}.sigma2(burn_in:end,1),'ko');
        subplot(2,2,4);
        plot(all_results{jj}{ii}.sigma2(burn_in:end,2),'ko');
        waitforbuttonpress;
    end
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
