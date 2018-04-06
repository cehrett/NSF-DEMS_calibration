function log_likelihood = resample_likelihood(desired_obs,theta,sigma2,...
    settings,Sigma_eta_xx)

% Get control settings for desired data
obs_x = settings.obs_x;

% Number of desired observations
n = size(obs_x,1);
num_out = length(desired_obs);
num_obs = n/num_out; % Total # of distinct multivariate obs

% Get nugget-making function
nugsize = settings.nugsize;

% Get simulation observations
eta = settings.eta;

% Get standardized desired observation values 
sim_output_means = mean(settings.output_means') ;
sim_output_sds   = mean(settings.output_sds')   ;
y = (desired_obs - sim_output_means) ./ sim_output_sds ;
y = repelem(y,size(obs_x,1)/length(y))' ;

% Make vector of desired and simulated observations
z = [ y ; eta ];

% Get prior on theta, and related functions/settings
log_theta_prior = settings.log_theta_prior;
Cost_lambda = settings.Cost_lambda;

% Get prior on sigma2
log_sigma2_prior = settings.log_sigma2_prior ;

% Get control and calibration settings at simulation observations
num_cntrl = size(obs_x,2) ;
sim_x = settings.sim_xt(:,1:num_cntrl) ;
sim_t = settings.sim_xt(:,num_cntrl+1:end) ;

% Omega, rho, lambda covariance parameters
omega=settings.omega; 
rho=settings.rho; 
lambda=settings.lambda;

% Arrange theta sample for use with n observations
obs_theta = repmat(theta,n,1) ; 

% Get Sigma_eta_yy (cov matrix for (desired) observations):
Sigma_eta_yy = gp_cov(omega,obs_x,obs_x,rho,...
    obs_theta,obs_theta,lambda,false);

% Get Sigma_eta_xy, and hence Sigma_eta_yx
Sigma_eta_xy = gp_cov(omega,sim_x,obs_x,...
    rho,sim_t,obs_theta,lambda,false);
Sigma_eta_yx = Sigma_eta_xy';

% Get Sigma_y
sigma2_long = repelem(sigma2,num_obs);
Sigma_y = diag(sigma2_long);

% Combine these to get Sigma_z
Sigma_z = [ Sigma_eta_yy + Sigma_y  Sigma_eta_yx ; Sigma_eta_xy ...
    Sigma_eta_xx ] ;

% Add a nugget to Sigma_z for computational stability
Sigma_z = Sigma_z + eye(size(Sigma_z)) * nugsize(Sigma_z);
% Get log likelihood of theta
L_D = logmvnpdf(z',0,Sigma_z);
log_likelihood = L_D + log_theta_prior(theta,Cost_lambda) + ...
    log_sigma2_prior(sigma2);

end