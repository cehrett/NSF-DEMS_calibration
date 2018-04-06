function [post_mean,weights] = resample_mean(settings, results , ...
    new_desired_obs)

% Get total number of samples
M=settings.M - settings.burn_in;

% Get old desired obs
old_desired_obs = results.desired_obs;

% Get samples of theta and sigma_2
theta = results.samples(settings.burn_in:end,:);
sigma2 = results.sigma2(settings.burn_in:end,:);

% This will store individual weights
weights = zeros(M,1);

% Get omega, rho, lambda
omega = settings.omega;
rho = settings.rho;
lambda = settings.lambda;

% Get control and calibration settings at simulation observations
num_cntrl = size(settings.obs_x,2) ;
sim_x = settings.sim_xt(:,1:num_cntrl) ;
sim_t = settings.sim_xt(:,num_cntrl+1:end) ;

% This is the covariance matrix for the simulation observation points.
% We'll need it as part of a larger cov matrix every time we get a
% likelihood, so we should calculate it once here and reuse it well.
Sigma_eta_xx = gp_cov(omega,sim_x,sim_x,rho,sim_t,sim_t,lambda,false);

% Call exterior function to get likelihood of original and new data at each
% of the M samples
msg=0; % For console output
for ii = 1 : M
    weight_num = resample_likelihood(new_desired_obs,theta(ii,:),...
        sigma2(ii,:),settings,Sigma_eta_xx);
    weight_den = resample_likelihood(old_desired_obs,theta(ii,:),...
        sigma2(ii,:),settings,Sigma_eta_xx);
    weights(ii) = weight_num/weight_den;
    
    % Console output to let us know progress
    if mod(ii,100) == 0
        fprintf(repmat('\b',1,msg));
        msg=fprintf('Completed: %d/%d',ii,M);
    end
end

% Calculate posterior mean
post_mean = sum(samples .* weights) / sum(weights);

end